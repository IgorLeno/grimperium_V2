import { useEffect, useState } from 'react'
import { postJson } from './api'
import type { TrainingFailure, TrainingInFlight, TrainingStatus } from './api'
import { PageHero, SectionHeading } from './StoryVisuals'
import { HEARTBEAT_ALERT_S, elapsedSince, formatDuration, heartbeatAge, runStats } from './trainingRun'
import { useApi } from './useApi'

const POLL_MS = 5000
const FEED_SIZE = 8
const number = new Intl.NumberFormat('pt-BR')
const percent = new Intl.NumberFormat('pt-BR', { maximumFractionDigits: 1 })
const kcal = new Intl.NumberFormat('pt-BR', { minimumFractionDigits: 1, maximumFractionDigits: 1 })
const rate = new Intl.NumberFormat('pt-BR', { maximumFractionDigits: 1 })

const JOB_LABEL: Record<string, string> = {
  idle: 'Parado',
  running: 'Em execução',
  stopping: 'Parando (terminando as moléculas em andamento)',
  finished: 'Concluído',
  failed: 'Falhou',
  stale: 'Interrompido sem aviso',
}

/** Relógio de 1 s para cronômetros e idade do heartbeat entre um poll e outro. */
function useNow(active: boolean): number {
  const [now, setNow] = useState(() => Date.now())
  useEffect(() => {
    if (!active) return undefined
    const id = window.setInterval(() => setNow(Date.now()), 1000)
    return () => window.clearInterval(id)
  }, [active])
  return now
}

export function TrainingView() {
  const { state, reload } = useApi<TrainingStatus>('/api/training/status', { pollMs: POLL_MS })
  const [workers, setWorkers] = useState('3')
  const [threads, setThreads] = useState('4')
  const [limit, setLimit] = useState('')
  const [busy, setBusy] = useState(false)
  const [actionError, setActionError] = useState<string | null>(null)

  const act = async (path: string, body?: unknown) => {
    setBusy(true)
    setActionError(null)
    try {
      await postJson<TrainingStatus>(path, body)
      reload()
    } catch (error) {
      setActionError(String(error instanceof Error ? error.message : error))
    } finally {
      setBusy(false)
    }
  }

  return (
    <div className="view">
      <PageHero
        variant="control" eyebrow="Calcular PM7 · Semi-Imperium" title="Calcule as moléculas" accent="que ensinam o modelo."
        lead="Esta rotina leva cada molécula do SMILES até a entalpia de formação PM7 com o protocolo canônico e guarda a geometria otimizada. O painel só aciona e acompanha: o cálculo roda no back-end."
      />
      {state.status === 'loading' && <p role="status" className="loading-state">Carregando andamento…</p>}
      {state.status === 'error' && (
        <p role="alert" className="error-state">{state.message} <button type="button" onClick={reload}>Tentar de novo</button></p>
      )}
      {state.status === 'ready' && <Ready data={state.data} busy={busy} workers={workers} threads={threads} limit={limit}
        onWorkers={setWorkers} onThreads={setThreads} onLimit={setLimit} actionError={actionError}
        onStart={() => act('/api/training/start', {
          workers: Number(workers), crest_threads: Number(threads), limit: limit.trim() === '' ? null : Number(limit),
        })}
        onStop={() => act('/api/training/stop')} />}
    </div>
  )
}

interface ReadyProps {
  data: TrainingStatus
  busy: boolean
  workers: string
  threads: string
  limit: string
  actionError: string | null
  onWorkers: (value: string) => void
  onThreads: (value: string) => void
  onLimit: (value: string) => void
  onStart: () => void
  onStop: () => void
}

function Ready({ data, busy, workers, threads, limit, actionError, onWorkers, onThreads, onLimit, onStart, onStop }: ReadyProps) {
  const { job } = data
  const active = job.state === 'running' || job.state === 'stopping'
  const now = useNow(active)
  const canStart = !active && !busy

  return (
    <>
      {job.started_at && <RunPanel data={data} active={active} now={now} />}
      {active && <WorkersPanel data={data} now={now} />}

      <section className="panel" aria-label="Controle do cálculo">
        <SectionHeading eyebrow="Controle" title="Iniciar ou parar" />
        <div className="training-controls">
          <label>Moléculas em paralelo
            <input type="number" min={1} max={8} value={workers} onChange={(e) => onWorkers(e.target.value)} disabled={active} />
          </label>
          <label>Threads por molécula
            <input type="number" min={1} max={64} value={threads} onChange={(e) => onThreads(e.target.value)} disabled={active} />
          </label>
          <label>Limite de moléculas (opcional)
            <input type="number" min={1} value={limit} onChange={(e) => onLimit(e.target.value)} disabled={active} />
          </label>
          <button type="button" className="primary-button" disabled={!canStart} onClick={onStart}>Iniciar cálculo</button>
          <button type="button" className="ghost-button" disabled={!active || job.state === 'stopping' || busy} onClick={onStop}>Parar</button>
        </div>
        {!active && <CpuBudget workers={Number(workers)} threads={Number(threads)} cpus={data.cpu_count} />}
        {active && <p role="status">Em execução com {job.workers} em paralelo: {number.format(job.done)} concluídas, {number.format(job.failed)} com falha nesta rodada.</p>}
        {job.state === 'stale' && <p role="status">O último cálculo foi interrompido sem aviso; as moléculas em andamento voltam para a fila no próximo início.</p>}
        {job.message && !active && <p role="status">{job.message}</p>}
        {actionError && <p role="alert" className="error-state">{actionError}</p>}
      </section>

      <Feed data={data} />
      <QueuePanel data={data} />
      <p className="train-note">{number.format(data.geometries)} geometrias finais guardadas no banco do treinamento.</p>
    </>
  )
}

/** Workers × threads contra os núcleos: passar disso deixa tudo mais lento. */
function CpuBudget({ workers, threads, cpus }: { workers: number; threads: number; cpus: number }) {
  if (!(workers > 0 && threads > 0)) return null
  const used = workers * threads
  const over = used > cpus
  return (
    <p className={`cpu-budget${over ? ' cpu-budget--over' : ''}`}>
      CREST: {workers} × {threads} = <strong>{used}</strong> de {cpus} núcleos
      {over && ' — acima dos núcleos disponíveis, as moléculas vão disputar CPU e ficar mais lentas.'}
    </p>
  )
}

function RunPanel({ data, active, now }: { data: TrainingStatus; active: boolean; now: number }) {
  const stats = runStats(data, now)
  if (!stats) return null
  const { job } = data
  const share = stats.target === 0 ? 0 : (100 * stats.processed) / stats.target
  const okShare = stats.target === 0 ? 0 : (100 * job.done) / stats.target
  const age = heartbeatAge(data, now)
  const silent = active && age !== null && age > HEARTBEAT_ALERT_S
  const label = `Rodada: ${number.format(stats.processed)} de ${number.format(stats.target)} moléculas processadas`

  return (
    <section className="panel run-panel" aria-label="Rodada atual">
      <SectionHeading eyebrow={active ? 'Rodada atual' : 'Última rodada'} title={`${number.format(stats.processed)} de ${number.format(stats.target)} moléculas`}
        badge={JOB_LABEL[job.state] ?? job.state}>
        {job.limit ? `Limite de ${number.format(job.limit)} nesta rodada.` : 'Sem limite: vai até esvaziar a fila.'}
      </SectionHeading>
      <div className={`run-bar${active ? ' run-bar--active' : ''}`} role="progressbar" aria-label={label}
        aria-valuemin={0} aria-valuemax={stats.target} aria-valuenow={stats.processed}>
        <span className="run-bar__fail" style={{ width: `${share}%` }} />
        <span className="run-bar__ok" style={{ width: `${okShare}%` }} />
        <b>{percent.format(share)}%</b>
      </div>
      <dl className="run-stats">
        <div><dt>Decorrido</dt><dd>{formatDuration(stats.elapsedS)}</dd></div>
        <div><dt>{active ? 'Termina em' : 'Restantes'}</dt>
          <dd>{active
            ? (stats.etaS === null ? 'estimando…' : `~${formatDuration(stats.etaS)}`)
            : number.format(stats.target - stats.processed)}</dd></div>
        <div><dt>Ritmo</dt><dd>{stats.perHour === null ? '—' : `${rate.format(stats.perHour)} mol/h`}</dd></div>
        <div><dt>Tempo médio</dt><dd>{data.activity.mean_elapsed_s === null ? '—' : formatDuration(data.activity.mean_elapsed_s)}</dd></div>
        <div><dt>Falhas</dt><dd>{stats.processed === 0 ? '—' : `${number.format(job.failed)} (${percent.format((100 * job.failed) / stats.processed)}%)`}</dd></div>
      </dl>
      {active && age !== null && (
        <p className={`heartbeat${silent ? ' heartbeat--silent' : ''}`} role="status">
          {silent
            ? `Sem sinal do coordenador há ${formatDuration(age)}: pode estar travado.`
            : `Coordenador ativo · último sinal há ${formatDuration(age)}`}
        </p>
      )}
    </section>
  )
}

function WorkersPanel({ data, now }: { data: TrainingStatus; now: number }) {
  const slots = Math.max(data.job.workers, data.activity.in_flight.length)
  const cards: (TrainingInFlight | null)[] = Array.from({ length: slots }, (_, i) => data.activity.in_flight[i] ?? null)
  return (
    <section className="panel" aria-label="Moléculas em cálculo">
      <SectionHeading eyebrow="Agora" title="Em cálculo" />
      <ul className="worker-grid">
        {cards.map((molecule, index) => (
          <WorkerCard key={molecule?.mol_id ?? `free-${index}`} index={index} molecule={molecule}
            mean={data.activity.mean_elapsed_s} now={now} />
        ))}
      </ul>
    </section>
  )
}

function WorkerCard({ index, molecule, mean, now }: { index: number; molecule: TrainingInFlight | null; mean: number | null; now: number }) {
  if (!molecule) {
    return (
      <li className="worker-card worker-card--free">
        <span className="worker-card__slot">Worker {index + 1}</span>
        <p>Aguardando a próxima molécula…</p>
      </li>
    )
  }
  const elapsed = elapsedSince(molecule.started_at, now)
  // Comparado à média da rodada; fica em 95% enquanto não termina para não prometer o fim.
  const share = mean && elapsed !== null ? Math.min(95, (100 * elapsed) / mean) : null
  return (
    <li className="worker-card">
      <span className="worker-card__slot"><i className="worker-card__spinner" aria-hidden="true" />Worker {index + 1}</span>
      <strong>{molecule.mol_id}</strong>
      <code title={molecule.smiles}>{molecule.smiles}</code>
      <span className="worker-card__meta">{molecule.nheavy} átomos pesados · {elapsed === null ? '—' : formatDuration(elapsed)}</span>
      <div className={`worker-card__bar${share === null ? ' worker-card__bar--unknown' : ''}`} aria-hidden="true">
        {share !== null && <span style={{ width: `${share}%` }} />}
      </div>
      {share !== null && elapsed !== null && mean !== null && elapsed > mean && (
        <span className="worker-card__meta">Acima do tempo médio ({formatDuration(mean)})</span>
      )}
    </li>
  )
}

type FeedItem =
  | { kind: 'ok'; at: string; mol_id: string; smiles: string; pm7: number; delta: number | null; elapsed: number | null; conformers: number | null }
  | { kind: 'fail'; at: string; mol_id: string; reason: string }

function feedItems(data: TrainingStatus): FeedItem[] {
  const ok: FeedItem[] = data.activity.recent.map((m) => ({
    kind: 'ok', at: m.finished_at, mol_id: m.mol_id, smiles: m.smiles, pm7: m.h298_pm7,
    delta: m.delta_kcal_mol, elapsed: m.elapsed_s, conformers: m.conformers_evaluated,
  }))
  const failed: FeedItem[] = data.failures.map((f: TrainingFailure) => ({ kind: 'fail', at: f.at, mol_id: f.mol_id, reason: f.reason }))
  return [...ok, ...failed].sort((a, b) => Date.parse(b.at) - Date.parse(a.at)).slice(0, FEED_SIZE)
}

function Feed({ data }: { data: TrainingStatus }) {
  const items = feedItems(data)
  if (items.length === 0) return null
  return (
    <section className="panel" aria-label="Últimas moléculas">
      <SectionHeading eyebrow="Histórico" title="Últimas moléculas processadas" />
      <ul className="feed">
        {items.map((item) => item.kind === 'ok' ? (
          <li key={`ok-${item.mol_id}-${item.at}`} className="feed__item feed__item--ok">
            <span className="feed__mark" aria-label="Concluída">✓</span>
            <strong>{item.mol_id}</strong>
            <code title={item.smiles}>{item.smiles}</code>
            <span>ΔHf PM7 {kcal.format(item.pm7)} kcal/mol{item.delta !== null && ` · Δ ${kcal.format(item.delta)}`}</span>
            <small>{item.elapsed !== null && formatDuration(item.elapsed)}{item.conformers !== null && ` · ${item.conformers} ${item.conformers === 1 ? 'confôrmero' : 'confôrmeros'}`}</small>
          </li>
        ) : (
          <li key={`fail-${item.mol_id}-${item.at}`} className="feed__item feed__item--fail">
            <span className="feed__mark" aria-label="Falhou">✗</span>
            <strong>{item.mol_id}</strong>
            <span>{item.reason}</span>
          </li>
        ))}
      </ul>
    </section>
  )
}

function QueuePanel({ data }: { data: TrainingStatus }) {
  const { progress } = data
  const queued = progress.pending + progress.rerun
  const segments = [
    { key: 'ok', label: 'Concluídas', count: progress.ok, className: 'training-bar__ok' },
    { key: 'run', label: 'Em execução', count: progress.running, className: 'training-bar__run' },
    { key: 'queue', label: 'Na fila', count: queued, className: 'training-bar__queue' },
    { key: 'skip', label: 'Descartadas', count: progress.skipped, className: 'training-bar__skip' },
  ]
  const share = (count: number) => (progress.total === 0 ? 0 : (100 * count) / progress.total)
  const summary = `Andamento do treinamento: ${segments.map((s) => `${s.label}: ${number.format(s.count)} de ${number.format(progress.total)}`).join('; ')}`

  return (
    <section className="panel" aria-label="Andamento do treinamento">
      <SectionHeading eyebrow="Fila completa" title={`${number.format(progress.ok)} de ${number.format(progress.total)} moléculas concluídas`}>
        Fila de treinamento: uma molécula por SMILES da referência CBS, todas no protocolo 2; {number.format(queued)} ainda na fila.
      </SectionHeading>
      <div className="training-bar training-bar--thin" role="img" aria-label={summary}>
        {segments.map((s) => s.count > 0 && <span key={s.key} className={s.className} style={{ width: `${share(s.count)}%` }} />)}
      </div>
      <ul className="training-legend">
        {segments.map((s) => (
          <li key={s.key}><i className={s.className} aria-hidden="true" />
            {s.key === 'queue' ? `${number.format(s.count)} na fila` : `${s.label}: ${number.format(s.count)}`}
            {' '}<small>({percent.format(share(s.count))}% de {number.format(progress.total)})</small>
          </li>
        ))}
      </ul>
    </section>
  )
}

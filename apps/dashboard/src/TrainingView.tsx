import { useState } from 'react'
import { postJson } from './api'
import type { TrainingStatus } from './api'
import { PageHero, SectionHeading } from './StoryVisuals'
import { useApi } from './useApi'

const POLL_MS = 5000
const number = new Intl.NumberFormat('pt-BR')
const percent = new Intl.NumberFormat('pt-BR', { maximumFractionDigits: 1 })

const JOB_LABEL: Record<string, string> = {
  idle: 'Parado',
  running: 'Em execução',
  stopping: 'Parando (terminando as moléculas em andamento)',
  finished: 'Concluído',
  failed: 'Falhou',
  stale: 'Interrompido sem aviso',
}

export function TrainingView() {
  const { state, reload } = useApi<TrainingStatus>('/api/training/status', { pollMs: POLL_MS })
  const [workers, setWorkers] = useState('3')
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
        variant="control" eyebrow="Treinamento · Semi-Imperium" title="Calcule as moléculas" accent="que ensinam o modelo."
        lead="Esta rotina leva cada molécula do SMILES até a entalpia de formação PM7 com o protocolo canônico e guarda a geometria otimizada. O painel só aciona e acompanha: o cálculo roda no back-end."
      />
      {state.status === 'loading' && <p role="status" className="loading-state">Carregando andamento…</p>}
      {state.status === 'error' && (
        <p role="alert" className="error-state">{state.message} <button type="button" onClick={reload}>Tentar de novo</button></p>
      )}
      {state.status === 'ready' && <Ready data={state.data} busy={busy} workers={workers} limit={limit}
        onWorkers={setWorkers} onLimit={setLimit} actionError={actionError}
        onStart={() => act('/api/training/start', { workers: Number(workers), limit: limit.trim() === '' ? null : Number(limit) })}
        onStop={() => act('/api/training/stop')} />}
    </div>
  )
}

interface ReadyProps {
  data: TrainingStatus
  busy: boolean
  workers: string
  limit: string
  actionError: string | null
  onWorkers: (value: string) => void
  onLimit: (value: string) => void
  onStart: () => void
  onStop: () => void
}

function Ready({ data, busy, workers, limit, actionError, onWorkers, onLimit, onStart, onStop }: ReadyProps) {
  const { progress, job } = data
  const queued = progress.pending + progress.rerun
  const segments = [
    { key: 'ok', label: 'Concluídas', count: progress.ok, className: 'training-bar__ok' },
    { key: 'run', label: 'Em execução', count: progress.running, className: 'training-bar__run' },
    { key: 'queue', label: 'Na fila', count: queued, className: 'training-bar__queue' },
    { key: 'skip', label: 'Descartadas', count: progress.skipped, className: 'training-bar__skip' },
  ]
  const active = job.state === 'running' || job.state === 'stopping'
  const share = (count: number) => (progress.total === 0 ? 0 : (100 * count) / progress.total)
  const summary = `Andamento do treinamento: ${segments.map((s) => `${s.label}: ${number.format(s.count)} de ${number.format(progress.total)}`).join('; ')}`
  const canStart = !active && !busy

  return (
    <>
      <section className="panel" aria-label="Andamento do treinamento">
        <SectionHeading eyebrow="Andamento" title={`${number.format(progress.ok)} de ${number.format(progress.total)} moléculas concluídas`}
          badge={JOB_LABEL[job.state] ?? job.state}>
          Fila de treinamento: uma molécula por SMILES da referência CBS, todas no protocolo 2; {number.format(queued)} ainda na fila.
        </SectionHeading>
        <div className="training-bar" role="img" aria-label={summary}>
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

      <section className="panel" aria-label="Controle do treinamento">
        <SectionHeading eyebrow="Controle" title="Iniciar ou parar" />
        <div className="training-controls">
          <label>Moléculas em paralelo
            <input type="number" min={1} max={8} value={workers} onChange={(e) => onWorkers(e.target.value)} disabled={active} />
          </label>
          <label>Limite de moléculas (opcional)
            <input type="number" min={1} value={limit} onChange={(e) => onLimit(e.target.value)} disabled={active} />
          </label>
          <button type="button" className="primary-button" disabled={!canStart} onClick={onStart}>Iniciar treinamento</button>
          <button type="button" className="ghost-button" disabled={!active || job.state === 'stopping' || busy} onClick={onStop}>Parar</button>
        </div>
        {active && <p role="status">Em execução com {job.workers} em paralelo: {number.format(job.done)} concluídas, {number.format(job.failed)} com falha nesta rodada.</p>}
        {job.state === 'stale' && <p role="status">O último treinamento foi interrompido sem aviso; as moléculas em andamento voltam para a fila no próximo início.</p>}
        {job.message && !active && <p role="status">{job.message}</p>}
        {actionError && <p role="alert" className="error-state">{actionError}</p>}
      </section>

      {data.failures.length > 0 && (
        <section className="panel" aria-label="Falhas recentes">
          <SectionHeading eyebrow="Falhas" title="Últimas moléculas que não passaram" />
          <ul className="training-failures">
            {data.failures.map((f) => (
              <li key={`${f.mol_id}-${f.at}`}><strong>{f.mol_id}</strong> <span>{f.reason}</span></li>
            ))}
          </ul>
        </section>
      )}
      <p className="train-note">{number.format(data.geometries)} geometrias finais guardadas no banco do treinamento.</p>
    </>
  )
}

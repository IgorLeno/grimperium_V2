import { useEffect, useState } from 'react'
import { downloadCellPackage, getJson, postJson } from './api'
import type { CellCreated, CellSummary, TrainingInFlight, TrainingStatus } from './api'
import { Callout, SectionHeading } from './StoryVisuals'
import { HEARTBEAT_ALERT_S, elapsedSince, formatDuration, heartbeatAge } from './trainingRun'
import { useApi } from './useApi'
import { useNow } from './useNow'

const POLL_MS = 10000
const CENTRAL = 'central'
// Comando do guia (docs/guia-celulas.md, seção 7): o dashboard não religa máquina nenhuma.
const RESTART_COMMAND = "wsl --shutdown; Start-ScheduledTask -TaskName 'Grimperium Celula'"

const STATE_LABEL: Record<string, string> = {
  online: 'Online',
  late: 'Atrasada',
  silent: 'Sem sinal',
  never_seen: 'Ainda não conectou',
  version_mismatch: 'Versão divergente',
  blocked: 'Bloqueada',
  revoked: 'Revogada',
}

const STATE_TONE: Record<string, string> = {
  online: 'good',
  late: 'warning',
  silent: 'critical',
  version_mismatch: 'critical',
  blocked: 'critical',
}

const JOB_LABEL: Record<string, string> = {
  idle: 'Parada',
  running: 'Calculando',
  stopping: 'Parando',
  finished: 'Parada',
  failed: 'Parada (última rodada falhou)',
  stale: 'Interrompida sem aviso',
}

/** O que ocupa um worker: molécula da fila, referência da admissão ou nada. */
type Slot = { kind: 'molecule'; molecule: TrainingInFlight } | { kind: 'admission' } | { kind: 'free' }

/**
 * Um bloco por computador: esta central e cada célula, com os workers, o
 * último sinal e as ações. Criar célula e baixar o pacote ficam no fim.
 */
export function CellsPanel({ training, onRefreshTraining }: { training: TrainingStatus; onRefreshTraining?: () => void }) {
  const { state, refresh } = useApi<CellSummary[]>('/api/cells', { pollMs: POLL_MS })
  const now = useNow(true)
  const [name, setName] = useState('')
  const [centralUrl, setCentralUrl] = useState('')
  const [busy, setBusy] = useState(false)
  const [message, setMessage] = useState<string | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    const controller = new AbortController()
    getJson<{ central_url: string }>('/api/cells/package-defaults', controller.signal).then(
      (data) => setCentralUrl((current) => current || (typeof data.central_url === 'string' ? data.central_url : '')),
      () => undefined, // sem sugestão: o campo fica vazio para digitar
    )
    return () => controller.abort()
  }, [])

  const run = async (action: () => Promise<string | null>): Promise<boolean> => {
    setBusy(true)
    setError(null)
    setMessage(null)
    try {
      setMessage(await action())
      refresh()
      onRefreshTraining?.()
      return true
    } catch (failure) {
      setError(String(failure instanceof Error ? failure.message : failure))
      return false
    } finally {
      setBusy(false)
    }
  }

  const url = centralUrl.trim()
  const path = (cell: CellSummary, action: string) => `/api/cells/${encodeURIComponent(cell.cell_id)}/${action}`
  const createAndDownload = () => run(async () => {
    const created = await postJson<CellCreated>('/api/cells', { name: name.trim() })
    const file = await downloadCellPackage(created.cell.cell_id, url)
    setName('')
    return `Célula "${created.cell.name}" criada; pacote ${file} baixado.`
  })
  const actions: CellActions = {
    busy,
    canPackage: url !== '',
    newPackage: (cell) => {
      if (!window.confirm(`Gerar um pacote novo para "${cell.name}"? O token muda: a instalação atual dessa célula para de funcionar até receber o pacote novo.`)) return
      void run(async () => `Pacote ${await downloadCellPackage(cell.cell_id, url)} baixado (token novo).`)
    },
    setAdmission: (cell, admission) => void run(async () => {
      await postJson<CellSummary>(path(cell, 'admission'), { admission })
      return null
    }),
    rename: (cell, newName) => run(async () => {
      await postJson<CellSummary>(path(cell, 'rename'), { name: newName })
      return `Célula renomeada para "${newName}".`
    }),
    release: (cell) => {
      if (!window.confirm(`Devolver à fila as ${cell.held} moléculas reservadas a "${cell.name}"? Outras células podem pegá-las. Se o resultado dela chegar depois, ainda é aceito enquanto ninguém concluir a molécula antes.`)) return
      void run(async () => {
        const { released } = await postJson<{ released: string[] }>(path(cell, 'release'))
        return `${released.length} ${released.length === 1 ? 'molécula devolvida' : 'moléculas devolvidas'} à fila.`
      })
    },
    revoke: (cell) => {
      if (!window.confirm(`Revogar "${cell.name}"? O token deixa de valer e a célula não recebe mais trabalho. Não dá para desfazer: para voltar, crie outra célula.`)) return
      void run(async () => {
        await postJson<CellSummary>(path(cell, 'revoke'))
        return `Célula "${cell.name}" revogada.`
      })
    },
    check: () => {
      refresh()
      onRefreshTraining?.()
    },
  }

  const cells = state.status === 'ready' && Array.isArray(state.data) ? state.data : []
  const inFlight = training.activity.in_flight
  const mean = training.activity.mean_elapsed_s

  return (
    <section className="panel" aria-label="Computadores">
      <SectionHeading eyebrow="Agora" title="Computadores">
        Esta central e as células: PCs com Windows + WSL2 que pedem moléculas a esta central pelo Tailscale. Guia: docs/guia-celulas.md.
      </SectionHeading>
      {message && <p role="status">{message}</p>}
      {error && <p role="alert" className="error-state">{error}</p>}
      {state.status === 'error' && <p role="alert" className="error-state">{state.message}</p>}
      <ul className="computers">
        <CentralRow training={training} inFlight={inFlight.filter((m) => !m.cell_id || m.cell_id === CENTRAL)} mean={mean} now={now} />
        {cells.map((cell) => (
          <CellRow key={cell.cell_id} cell={cell} inFlight={inFlight.filter((m) => m.cell_id === cell.cell_id)}
            mean={mean} now={now} actions={actions} />
        ))}
      </ul>

      <h3 className="computers__subtitle">Nova célula</h3>
      <Callout label="Atenção" tone="warning">
        O pacote da célula contém o token de acesso dela (config.json). Não envie por e-mail nem deixe em pasta
        compartilhada; apague o zip depois de instalar. Cada pacote gera um token novo e invalida o anterior.
      </Callout>
      <div className="training-controls">
        <label>Nome do computador
          <input value={name} onChange={(e) => setName(e.target.value)} placeholder="PC do laboratório" maxLength={80} />
        </label>
        <label>Endereço da central (Tailscale)
          <input value={centralUrl} onChange={(e) => setCentralUrl(e.target.value)} placeholder="http://100.x.y.z:8100" />
        </label>
        <button type="button" className="primary-button" disabled={busy || name.trim() === '' || url === ''}
          onClick={() => void createAndDownload()}>Criar célula e baixar pacote</button>
      </div>
    </section>
  )
}

interface CellActions {
  busy: boolean
  canPackage: boolean
  newPackage: (cell: CellSummary) => void
  setAdmission: (cell: CellSummary, admission: string) => void
  rename: (cell: CellSummary, name: string) => Promise<boolean>
  release: (cell: CellSummary) => void
  revoke: (cell: CellSummary) => void
  check: () => void
}

function CentralRow({ training, inFlight, mean, now }: { training: TrainingStatus; inFlight: TrainingInFlight[]; mean: number | null; now: number }) {
  const { job } = training
  const active = job.state === 'running' || job.state === 'stopping'
  const age = active ? heartbeatAge(training, now) : null
  const silent = age !== null && age > HEARTBEAT_ALERT_S
  const workers = active ? job.workers : 0
  const slots = slotsFor(inFlight, Math.max(workers, inFlight.length), 0)
  return (
    <li className="computer" aria-label="Esta central">
      <div className="computer__head">
        <strong className="computer__name">Esta central</strong>
        <span className={`computer__state computer__state--${active ? (silent ? 'critical' : 'good') : 'neutral'}`}>
          {JOB_LABEL[job.state] ?? job.state}
        </span>
        <small>{age === null ? 'Rotina local: controle em "Iniciar ou parar"' : `Último sinal há ${formatDuration(age)}`}</small>
        <small>{busyFree(inFlight.length, workers)}</small>
      </div>
      {slots.length > 0 && <WorkerGrid slots={slots} mean={mean} now={now} />}
    </li>
  )
}

function CellRow({ cell, inFlight, mean, now, actions }: { cell: CellSummary; inFlight: TrainingInFlight[]; mean: number | null; now: number; actions: CellActions }) {
  const [editing, setEditing] = useState(false)
  const [draft, setDraft] = useState(cell.name)
  const [checked, setChecked] = useState(false)
  const revoked = cell.state === 'revoked'
  // Só célula em contato recebe trabalho: nas outras, worker livre não está esperando nada.
  const working = cell.state === 'online' || cell.state === 'late'
  const workers = cell.workers ?? 0
  // `held` conta também as referências da admissão, que não estão na fila.
  const admissionRefs = Math.max(0, cell.held - inFlight.length)
  const slots = revoked ? [] : slotsFor(inFlight, working ? Math.max(workers, cell.held) : cell.held, admissionRefs)
  const age = elapsedSince(cell.last_seen_at, now)
  const save = async () => {
    const value = draft.trim()
    if (value === '' || value === cell.name) {
      setEditing(false)
      return
    }
    // Se a API recusar, o campo continua aberto com o que foi digitado.
    if (await actions.rename(cell, value)) setEditing(false)
  }

  return (
    <li className={`computer${revoked ? ' computer--revoked' : ''}`} aria-label={`Célula ${cell.name}`}>
      <div className="computer__head">
        {editing ? (
          <form className="computer__rename" onSubmit={(e) => { e.preventDefault(); void save() }}>
            <input aria-label="Novo nome" value={draft} maxLength={80} autoFocus onChange={(e) => setDraft(e.target.value)} />
            <button type="submit" className="ghost-button ghost-button--small" disabled={actions.busy || draft.trim() === ''}>Salvar</button>
            <button type="button" className="ghost-button ghost-button--small" onClick={() => { setDraft(cell.name); setEditing(false) }}>Cancelar</button>
          </form>
        ) : (
          <strong className="computer__name">{cell.name}</strong>
        )}
        <span className={`computer__state computer__state--${STATE_TONE[cell.state] ?? 'neutral'}`}>{STATE_LABEL[cell.state] ?? cell.state}</span>
        {cell.admission === 'pending' && !revoked && <span className="computer__state computer__state--info">Em admissão</span>}
        <small>{age === null ? 'Nenhum sinal ainda' : `Último sinal há ${formatDuration(age)}`}</small>
        {cell.host && <small>{cell.host} · {cell.workers ?? '?'} × {cell.crest_threads ?? '?'} threads</small>}
        {working && <small>{busyFree(cell.held, workers)}</small>}
        {!working && !revoked && cell.held > 0 && <small>{cell.held} {cell.held === 1 ? 'molécula reservada' : 'moléculas reservadas'}</small>}
      </div>

      {cell.state === 'version_mismatch' && (
        <p className="computer__alert" role="note">Versão divergente: a central recusou o registro. {cell.refusal ?? ''} Baixe um pacote novo e reinstale.</p>
      )}
      {cell.state !== 'version_mismatch' && cell.refusal && <p className="computer__note">Recusada: {cell.refusal}</p>}
      {cell.admission_note && (
        <p className={cell.admission === 'blocked' ? 'computer__alert' : 'computer__note'}>
          {cell.admission === 'blocked' ? `Bloqueada pelo teste de admissão: ${cell.admission_note}` : cell.admission_note}
        </p>
      )}

      {slots.length > 0 && <WorkerGrid slots={slots} mean={mean} now={now} />}

      {checked && !revoked && (
        <p className="computer__note">
          {cell.state === 'online'
            ? `Respondendo: último sinal há ${age === null ? '—' : formatDuration(age)}.`
            : 'Sem resposta recente. Pelo TeamViewer, no PowerShell da célula:'}
          {cell.state !== 'online' && <> <code>{RESTART_COMMAND}</code> (ou reinicie o Windows).</>}
        </p>
      )}

      {!revoked && (
        <div className="computer__actions">
          <button type="button" className="ghost-button ghost-button--small" onClick={() => { setChecked(true); actions.check() }}>Verificar agora</button>
          <button type="button" className="ghost-button ghost-button--small" disabled={actions.busy || cell.held === 0}
            onClick={() => actions.release(cell)}>Devolver à fila</button>
          {!editing && (
            <button type="button" className="ghost-button ghost-button--small" disabled={actions.busy}
              onClick={() => { setDraft(cell.name); setEditing(true) }}>Renomear</button>
          )}
          <button type="button" className="ghost-button ghost-button--small" disabled={actions.busy || !actions.canPackage}
            onClick={() => actions.newPackage(cell)}>Baixar pacote novo</button>
          {cell.admission !== 'passed' && (
            <button type="button" className="ghost-button ghost-button--small" disabled={actions.busy}
              onClick={() => actions.setAdmission(cell, 'passed')}>Admitir</button>
          )}
          {cell.admission !== 'blocked' && (
            <button type="button" className="ghost-button ghost-button--small" disabled={actions.busy}
              onClick={() => actions.setAdmission(cell, 'blocked')}>Bloquear</button>
          )}
          {cell.admission !== 'pending' && (
            <button type="button" className="ghost-button ghost-button--small" disabled={actions.busy}
              onClick={() => actions.setAdmission(cell, 'pending')}>Refazer admissão</button>
          )}
          <button type="button" className="ghost-button ghost-button--small ghost-button--danger" disabled={actions.busy}
            onClick={() => actions.revoke(cell)}>Revogar</button>
        </div>
      )}
    </li>
  )
}

function busyFree(busy: number, workers: number): string {
  if (workers === 0) return busy === 0 ? 'Nenhum worker' : `${busy} em cálculo`
  return `${busy} de ${workers} workers ocupados · ${Math.max(0, workers - busy)} livres`
}

/** Moléculas da fila, depois referências da admissão, depois workers livres. */
function slotsFor(inFlight: TrainingInFlight[], total: number, admissionRefs: number): Slot[] {
  const slots: Slot[] = inFlight.map((molecule) => ({ kind: 'molecule', molecule }))
  for (let i = 0; i < admissionRefs; i++) slots.push({ kind: 'admission' })
  while (slots.length < total) slots.push({ kind: 'free' })
  return slots
}

function WorkerGrid({ slots, mean, now }: { slots: Slot[]; mean: number | null; now: number }) {
  return (
    <ul className="worker-grid">
      {slots.map((slot, index) => (
        <WorkerCard key={slot.kind === 'molecule' ? slot.molecule.mol_id : `${slot.kind}-${index}`} index={index} slot={slot} mean={mean} now={now} />
      ))}
    </ul>
  )
}

function WorkerCard({ index, slot, mean, now }: { index: number; slot: Slot; mean: number | null; now: number }) {
  if (slot.kind === 'free') {
    return (
      <li className="worker-card worker-card--free">
        <span className="worker-card__slot">Worker {index + 1}</span>
        <p>Aguardando a próxima molécula…</p>
      </li>
    )
  }
  if (slot.kind === 'admission') {
    return (
      <li className="worker-card">
        <span className="worker-card__slot"><i className="worker-card__spinner" aria-hidden="true" />Worker {index + 1}</span>
        <strong>Referência de admissão</strong>
        <span className="worker-card__meta">Molécula já aceita na central, recalculada para comparar o HOF</span>
      </li>
    )
  }
  const { molecule } = slot
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

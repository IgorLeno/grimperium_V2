import { useEffect, useState } from 'react'
import { downloadCellPackage, getJson, postJson } from './api'
import type { CellCreated, CellSummary } from './api'
import { Callout, SectionHeading } from './StoryVisuals'
import { useApi } from './useApi'

const POLL_MS = 10000

const STATE_LABEL: Record<string, string> = {
  online: 'Online',
  late: 'Atrasada',
  silent: 'Sem sinal',
  never_seen: 'Ainda não conectou',
  version_mismatch: 'Versão divergente',
  blocked: 'Bloqueada',
  revoked: 'Revogada',
}

const ADMISSION_LABEL: Record<string, string> = {
  pending: 'em admissão',
  passed: 'admitida',
  blocked: 'bloqueada',
}

/**
 * Computadores remotos (células): criar, baixar o pacote de instalação e admitir.
 * A grade completa por computador fica para o Lote 5; aqui só o necessário
 * para instalar uma célula nova.
 */
export function CellsPanel() {
  const { state, reload } = useApi<CellSummary[]>('/api/cells', { pollMs: POLL_MS })
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

  const run = async (action: () => Promise<string | null>) => {
    setBusy(true)
    setError(null)
    setMessage(null)
    try {
      setMessage(await action())
      reload()
    } catch (failure) {
      setError(String(failure instanceof Error ? failure.message : failure))
    } finally {
      setBusy(false)
    }
  }

  const url = centralUrl.trim()
  const createAndDownload = () => run(async () => {
    const created = await postJson<CellCreated>('/api/cells', { name: name.trim() })
    const file = await downloadCellPackage(created.cell.cell_id, url)
    setName('')
    return `Célula "${created.cell.name}" criada; pacote ${file} baixado.`
  })
  const newPackage = (cell: CellSummary) => {
    if (!window.confirm(`Gerar um pacote novo para "${cell.name}"? O token muda: a instalação atual dessa célula para de funcionar até receber o pacote novo.`)) return
    void run(async () => `Pacote ${await downloadCellPackage(cell.cell_id, url)} baixado (token novo).`)
  }
  const setAdmission = (cell: CellSummary, admission: string) => run(async () => {
    await postJson<CellSummary>(`/api/cells/${encodeURIComponent(cell.cell_id)}/admission`, { admission })
    return null
  })

  const cells = state.status === 'ready' && Array.isArray(state.data) ? state.data : []

  return (
    <section className="panel" aria-label="Células">
      <SectionHeading eyebrow="Outros computadores" title="Células">
        Cada célula é um PC com Windows + WSL2 que pede moléculas a esta central pelo Tailscale. Guia: docs/guia-celulas.md.
      </SectionHeading>
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
      {message && <p role="status">{message}</p>}
      {error && <p role="alert" className="error-state">{error}</p>}
      {state.status === 'error' && <p role="alert" className="error-state">{state.message}</p>}
      {cells.length > 0 && (
        <ul className="training-failures cells-list">
          {cells.map((cell) => (
            <li key={cell.cell_id}>
              <strong>{cell.name}</strong>
              <span>{STATE_LABEL[cell.state] ?? cell.state} · {ADMISSION_LABEL[cell.admission] ?? cell.admission}</span>
              {cell.host && <small>{cell.host} · {cell.workers ?? '?'} × {cell.crest_threads ?? '?'} threads · {cell.held} em cálculo</small>}
              {cell.refusal && <small>Recusada: {cell.refusal}</small>}
              {cell.admission_note && (
                <small>{cell.admission === 'blocked' ? `Bloqueada pelo teste de admissão: ${cell.admission_note}` : cell.admission_note}</small>
              )}
              {cell.state !== 'revoked' && (
                <>
                  <button type="button" className="ghost-button ghost-button--small" disabled={busy || url === ''}
                    onClick={() => newPackage(cell)}>Baixar pacote novo</button>
                  {cell.admission !== 'passed' && (
                    <button type="button" className="ghost-button ghost-button--small" disabled={busy}
                      onClick={() => void setAdmission(cell, 'passed')}>Admitir</button>
                  )}
                  {cell.admission !== 'blocked' && (
                    <button type="button" className="ghost-button ghost-button--small" disabled={busy}
                      onClick={() => void setAdmission(cell, 'blocked')}>Bloquear</button>
                  )}
                  {cell.admission !== 'pending' && (
                    <button type="button" className="ghost-button ghost-button--small" disabled={busy}
                      onClick={() => void setAdmission(cell, 'pending')}>Refazer admissão</button>
                  )}
                </>
              )}
            </li>
          ))}
        </ul>
      )}
    </section>
  )
}

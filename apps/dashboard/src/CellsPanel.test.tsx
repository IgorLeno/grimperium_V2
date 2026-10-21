import { cleanup, fireEvent, render, screen, waitFor, within } from '@testing-library/react'
import { afterEach, expect, test, vi } from 'vitest'
import { CellsPanel } from './CellsPanel'
import type { CellSummary, TrainingStatus } from './api'

const cell: CellSummary = {
  cell_id: 'cell-abc', name: 'PC do Lab', state: 'never_seen', admission: 'pending', admission_note: null, host: null, cores: null,
  workers: null, crest_threads: null, program_versions: {}, refusal: null, created_at: '2026-10-21T00:00:00Z',
  registered_at: null, last_seen_at: null, revoked_at: null, held: 0,
}

const idle: TrainingStatus['job'] = { state: 'idle', pid: null, started_at: null, heartbeat_at: null, finished_at: null, workers: 0, limit: null, done: 0, failed: 0, message: '', stop_requested: false }

function training(overrides: Partial<TrainingStatus['activity']> = {}, job: TrainingStatus['job'] = idle): TrainingStatus {
  return {
    progress: { total: 10, ok: 0, pending: 10, rerun: 0, running: 0, skipped: 0 }, job, geometries: 0, failures: [],
    activity: { in_flight: [], recent: [], mean_elapsed_s: null, ...overrides }, cpu_count: 16,
  }
}

const ago = (s: number) => new Date(Date.now() - s * 1000).toISOString()

const online: CellSummary = {
  ...cell, state: 'online', admission: 'passed', host: 'lab-pc', cores: 8, workers: 2, crest_threads: 4,
  registered_at: ago(600), last_seen_at: ago(12), held: 1,
}

function json(data: unknown, code = 200): Response {
  return new Response(JSON.stringify(data), { status: code, headers: { 'Content-Type': 'application/json' } })
}

function mock(cells: CellSummary[]) {
  const calls: { path: string; method: string; body: string | null }[] = []
  vi.stubGlobal('fetch', vi.fn((input: RequestInfo | URL, init?: RequestInit) => {
    const path = String(input)
    const method = init?.method ?? 'GET'
    calls.push({ path, method, body: typeof init?.body === 'string' ? init.body : null })
    if (path === '/api/cells/package-defaults') return Promise.resolve(json({ central_url: 'http://100.64.0.1:8100' }))
    if (path === '/api/cells' && method === 'POST') return Promise.resolve(json({ cell, token: 'segredo' }))
    if (path.endsWith('/package')) {
      return Promise.resolve(new Response('zip', {
        status: 200, headers: { 'Content-Type': 'application/zip', 'Content-Disposition': 'attachment; filename="celula-pc-do-lab-abc.zip"' },
      }))
    }
    if (path.endsWith('/admission')) return Promise.resolve(json({ ...cell, admission: 'passed' }))
    if (path.endsWith('/rename')) return Promise.resolve(json({ ...cell, name: 'Novo nome' }))
    if (path.endsWith('/release')) return Promise.resolve(json({ released: ['mol_00007'] }))
    if (path.endsWith('/revoke')) return Promise.resolve(json({ ...cell, state: 'revoked' }))
    return Promise.resolve(json(cells))
  }))
  return calls
}

afterEach(() => { cleanup(); vi.unstubAllGlobals(); vi.restoreAllMocks() })

test('avisa que o pacote contém o token e sugere o endereço da central', async () => {
  mock([])
  render(<CellsPanel training={training()} />)
  expect(screen.getByText(/contém o token de acesso/)).toBeInTheDocument()
  await waitFor(() => expect(screen.getByLabelText('Endereço da central (Tailscale)')).toHaveValue('http://100.64.0.1:8100'))
})

test('cria a célula e baixa o pacote com o endereço informado', async () => {
  const calls = mock([])
  URL.createObjectURL = vi.fn(() => 'blob:x')
  URL.revokeObjectURL = vi.fn()
  const click = vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(() => undefined)
  render(<CellsPanel training={training()} />)
  await waitFor(() => expect(screen.getByLabelText('Endereço da central (Tailscale)')).toHaveValue('http://100.64.0.1:8100'))
  fireEvent.change(screen.getByLabelText('Nome do computador'), { target: { value: 'PC do Lab' } })
  fireEvent.click(screen.getByRole('button', { name: 'Criar célula e baixar pacote' }))
  expect(await screen.findByText(/pacote celula-pc-do-lab-abc\.zip baixado/)).toBeInTheDocument()
  const create = calls.find((c) => c.path === '/api/cells' && c.method === 'POST')
  expect(create?.body).toBe(JSON.stringify({ name: 'PC do Lab' }))
  const pkg = calls.find((c) => c.path === '/api/cells/cell-abc/package')
  expect(pkg?.body).toBe(JSON.stringify({ central_url: 'http://100.64.0.1:8100' }))
  expect(click).toHaveBeenCalled()
  // O token da criação nunca aparece na tela.
  expect(screen.queryByText(/segredo/)).not.toBeInTheDocument()
})

test('pacote novo pede confirmação e admitir chama a rota de admissão', async () => {
  const calls = mock([cell])
  const confirm = vi.spyOn(window, 'confirm').mockReturnValue(false)
  render(<CellsPanel training={training()} />)
  await screen.findByText('PC do Lab')
  await waitFor(() => expect(screen.getByLabelText('Endereço da central (Tailscale)')).toHaveValue('http://100.64.0.1:8100'))
  fireEvent.click(screen.getByRole('button', { name: 'Baixar pacote novo' }))
  expect(confirm).toHaveBeenCalledWith(expect.stringMatching(/token muda/))
  expect(calls.some((c) => c.path.endsWith('/package'))).toBe(false)
  fireEvent.click(screen.getByRole('button', { name: 'Admitir' }))
  await waitFor(() => expect(calls.find((c) => c.path === '/api/cells/cell-abc/admission')?.body).toBe(JSON.stringify({ admission: 'passed' })))
})

test('célula bloqueada pelo teste mostra o motivo e refazer volta para pending', async () => {
  const blocked: CellSummary = {
    ...cell, state: 'blocked', admission: 'blocked',
    admission_note: 'cbs_09698 (CO): HOF -48.70000 contra -48.93880 da central (Δ 0.23880 > 0.1 kcal/mol)',
  }
  const calls = mock([blocked])
  render(<CellsPanel training={training()} />)
  expect(await screen.findByText(/Bloqueada pelo teste de admissão: cbs_09698 \(CO\)/)).toBeInTheDocument()
  expect(screen.queryByRole('button', { name: 'Bloquear' })).not.toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: 'Refazer admissão' }))
  await waitFor(() => expect(calls.find((c) => c.path === '/api/cells/cell-abc/admission')?.body).toBe(JSON.stringify({ admission: 'pending' })))
})

test('uma linha por computador: central e célula, cada uma com suas moléculas', async () => {
  mock([online])
  render(<CellsPanel training={training({
    in_flight: [
      { mol_id: 'mol_central', smiles: 'C', nheavy: 1, started_at: ago(30), cell_id: 'central' },
      { mol_id: 'mol_00007', smiles: 'CCO', nheavy: 3, started_at: ago(90), cell_id: 'cell-abc' },
    ],
  }, { ...idle, state: 'running', workers: 1, started_at: ago(600), heartbeat_at: ago(5) })} />)
  const row = await screen.findByRole('listitem', { name: 'Célula PC do Lab' })
  expect(within(row).getByText('mol_00007')).toBeInTheDocument()
  expect(within(row).getByText(/1 min 30 s/)).toBeInTheDocument()
  expect(within(row).getByText('Online')).toBeInTheDocument()
  expect(within(row).getByText(/Último sinal há 1\d s/)).toBeInTheDocument()
  expect(within(row).getByText('1 de 2 workers ocupados · 1 livres')).toBeInTheDocument()
  expect(within(row).getByText(/Aguardando a próxima molécula/)).toBeInTheDocument()
  expect(within(row).queryByText('mol_central')).not.toBeInTheDocument()
  const central = screen.getByRole('listitem', { name: 'Esta central' })
  expect(within(central).getByText('mol_central')).toBeInTheDocument()
  expect(within(central).getByText('Calculando')).toBeInTheDocument()
})

test('renomear edita o nome na própria linha', async () => {
  const calls = mock([online])
  render(<CellsPanel training={training()} />)
  const row = await screen.findByRole('listitem', { name: 'Célula PC do Lab' })
  fireEvent.click(within(row).getByRole('button', { name: 'Renomear' }))
  fireEvent.change(within(row).getByLabelText('Novo nome'), { target: { value: 'Novo nome' } })
  fireEvent.click(within(row).getByRole('button', { name: 'Salvar' }))
  await waitFor(() => expect(calls.find((c) => c.path === '/api/cells/cell-abc/rename')?.body).toBe(JSON.stringify({ name: 'Novo nome' })))
  expect(await screen.findByText(/renomeada para "Novo nome"/)).toBeInTheDocument()
})

test('devolver à fila e revogar pedem confirmação antes de chamar a API', async () => {
  const calls = mock([online])
  const confirm = vi.spyOn(window, 'confirm').mockReturnValueOnce(false).mockReturnValue(true)
  render(<CellsPanel training={training()} />)
  const row = await screen.findByRole('listitem', { name: 'Célula PC do Lab' })
  fireEvent.click(within(row).getByRole('button', { name: 'Devolver à fila' }))
  expect(confirm).toHaveBeenLastCalledWith(expect.stringMatching(/Devolver à fila as 1 moléculas/))
  expect(calls.some((c) => c.path.endsWith('/release'))).toBe(false)
  fireEvent.click(within(row).getByRole('button', { name: 'Devolver à fila' }))
  expect(await screen.findByText('1 molécula devolvida à fila.')).toBeInTheDocument()
  fireEvent.click(within(row).getByRole('button', { name: 'Revogar' }))
  expect(confirm).toHaveBeenLastCalledWith(expect.stringMatching(/Não dá para desfazer/))
  await waitFor(() => expect(calls.some((c) => c.path === '/api/cells/cell-abc/revoke' && c.method === 'POST')).toBe(true))
})

test('devolver fica desabilitado sem moléculas reservadas', async () => {
  mock([{ ...online, held: 0 }])
  render(<CellsPanel training={training()} />)
  const row = await screen.findByRole('listitem', { name: 'Célula PC do Lab' })
  expect(within(row).getByRole('button', { name: 'Devolver à fila' })).toBeDisabled()
})

test('versão divergente e admissão em andamento aparecem na linha', async () => {
  mock([
    { ...online, cell_id: 'cell-x', name: 'PC velho', state: 'version_mismatch', refusal: 'xtb 6.6.1 ≠ 6.7.1', held: 0 },
    { ...online, cell_id: 'cell-y', name: 'PC novo', admission: 'pending', admission_note: '2 de 5 referências conferidas', held: 2 },
  ])
  render(<CellsPanel training={training()} />)
  const old = await screen.findByRole('listitem', { name: 'Célula PC velho' })
  expect(within(old).getByRole('note')).toHaveTextContent(/Versão divergente.*xtb 6\.6\.1 ≠ 6\.7\.1/)
  // Recusada não recebe trabalho: nenhum worker "aguardando".
  expect(within(old).queryByText(/Aguardando a próxima molécula/)).not.toBeInTheDocument()
  const fresh = screen.getByRole('listitem', { name: 'Célula PC novo' })
  expect(within(fresh).getByText('Em admissão')).toBeInTheDocument()
  expect(within(fresh).getByText('2 de 5 referências conferidas')).toBeInTheDocument()
  // As referências não estão na fila: os dois workers mostram a admissão.
  expect(within(fresh).getAllByText('Referência de admissão')).toHaveLength(2)
})

test('verificar agora busca de novo e mostra o comando de reinício para célula sem sinal', async () => {
  const calls = mock([{ ...online, state: 'silent', last_seen_at: ago(900) }])
  const refreshTraining = vi.fn()
  render(<CellsPanel training={training()} onRefreshTraining={refreshTraining} />)
  const row = await screen.findByRole('listitem', { name: 'Célula PC do Lab' })
  expect(within(row).getByText('1 molécula reservada')).toBeInTheDocument()
  expect(within(row).queryByText(/Aguardando a próxima molécula/)).not.toBeInTheDocument()
  const before = calls.filter((c) => c.path === '/api/cells' && c.method === 'GET').length
  fireEvent.click(within(row).getByRole('button', { name: 'Verificar agora' }))
  expect(within(row).getByText(/Start-ScheduledTask -TaskName 'Grimperium Celula'/)).toBeInTheDocument()
  expect(refreshTraining).toHaveBeenCalled()
  await waitFor(() => expect(calls.filter((c) => c.path === '/api/cells' && c.method === 'GET').length).toBeGreaterThan(before))
})

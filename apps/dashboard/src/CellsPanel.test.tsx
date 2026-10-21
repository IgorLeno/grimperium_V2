import { cleanup, fireEvent, render, screen, waitFor } from '@testing-library/react'
import { afterEach, expect, test, vi } from 'vitest'
import { CellsPanel } from './CellsPanel'
import type { CellSummary } from './api'

const cell: CellSummary = {
  cell_id: 'cell-abc', name: 'PC do Lab', state: 'never_seen', admission: 'pending', admission_note: null, host: null, cores: null,
  workers: null, crest_threads: null, program_versions: {}, refusal: null, created_at: '2026-10-21T00:00:00Z',
  registered_at: null, last_seen_at: null, revoked_at: null, held: 0,
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
    return Promise.resolve(json(cells))
  }))
  return calls
}

afterEach(() => { cleanup(); vi.unstubAllGlobals(); vi.restoreAllMocks() })

test('avisa que o pacote contém o token e sugere o endereço da central', async () => {
  mock([])
  render(<CellsPanel />)
  expect(screen.getByText(/contém o token de acesso/)).toBeInTheDocument()
  await waitFor(() => expect(screen.getByLabelText('Endereço da central (Tailscale)')).toHaveValue('http://100.64.0.1:8100'))
})

test('cria a célula e baixa o pacote com o endereço informado', async () => {
  const calls = mock([])
  URL.createObjectURL = vi.fn(() => 'blob:x')
  URL.revokeObjectURL = vi.fn()
  const click = vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(() => undefined)
  render(<CellsPanel />)
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
  render(<CellsPanel />)
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
  render(<CellsPanel />)
  expect(await screen.findByText(/Bloqueada pelo teste de admissão: cbs_09698 \(CO\)/)).toBeInTheDocument()
  expect(screen.queryByRole('button', { name: 'Bloquear' })).not.toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: 'Refazer admissão' }))
  await waitFor(() => expect(calls.find((c) => c.path === '/api/cells/cell-abc/admission')?.body).toBe(JSON.stringify({ admission: 'pending' })))
})

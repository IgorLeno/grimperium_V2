import { cleanup, fireEvent, render, screen, within } from '@testing-library/react'
import type { ReactNode } from 'react'
import { afterEach, expect, test, vi } from 'vitest'
import { App } from './App'

vi.mock('recharts', () => {
  const container = ({ children }: { children?: ReactNode }) => <div>{children}</div>
  const empty = () => null
  return { CartesianGrid: empty, Legend: empty, ReferenceLine: empty,
    ResponsiveContainer: container, Scatter: empty, ScatterChart: container,
    Tooltip: empty, XAxis: empty, YAxis: empty }
})

const progress = {
  legacy: { protocol: 'crest_pm7', path: '/data/pm7.csv', total: 4,
    status: { OK: 2, Pending: 1, Failed: 1 }, cbs_quality_flag: { OK: 2, SUSPECT: 1, unset: 1 } },
  semi_imperium: { path: '/data/store', exists: true, molecules: 1, runs: 2,
    calculations: 3, states: { pending: 1, failed: 2 } },
}
const summary = { path: '/data/exp.csv', rows: 4,
  by_source: { atct: 2, 'nist:gas': 2 }, by_source_family: { atct: 2, nist: 2 },
  flags: { conflict_flag: 1, validation_only: 1, low_priority_only: 0 },
  uncertainty_kind: { reported: 3, none: 1 } }
const comparison = { path: '/data/vs_cbs.csv', units: 'kcal/mol',
  stats: { n: 2, mae: 2, bias: 2, rmse: 2.24 }, points: [
    { exp_id: 'exp_1', mol_id: 'cbs_1', smiles: 'C', nheavy: 1, source: 'atct', H298_exp: -17, H298_cbs: -16, cbs_minus_exp: 1, uncertainty: .1, conflict_flag: false, low_priority_only: false, validation_only: false },
    { exp_id: 'exp_2', mol_id: 'cbs_2', smiles: 'CCO', nheavy: 3, source: 'nist:gas', H298_exp: -56, H298_cbs: -53, cbs_minus_exp: 3, uncertainty: null, conflict_flag: true, low_priority_only: false, validation_only: true },
  ] }

function json(data: unknown, status = 200): Response {
  return new Response(JSON.stringify(data), { status, headers: { 'Content-Type': 'application/json' } })
}

function mockApi(overrides: Record<string, Response | Response[]> = {}) {
  const data: Record<string, Response | Response[]> = {
    '/api/health': json({ status: 'ok', read_only: true }),
    '/api/dataset/progress': json(progress),
    '/api/experimental/summary': json(summary),
    '/api/experimental/comparison': json(comparison),
    ...overrides,
  }
  const fetchMock = vi.fn((input: RequestInfo | URL) => {
    const reply = data[String(input)]
    return Promise.resolve(Array.isArray(reply) ? reply.shift() : reply)
  })
  vi.stubGlobal('fetch', fetchMock)
  return fetchMock
}

afterEach(() => { cleanup(); vi.unstubAllGlobals() })

test('mostra progresso real, protocolo legado e resumo do store', async () => {
  mockApi()
  render(<App />)
  expect(screen.getByText('Carregando progresso…')).toBeInTheDocument()
  expect(await screen.findByText('Protocolo legado: crest_pm7')).toBeInTheDocument()
  expect(screen.getByText('API ativa · somente leitura')).toBeInTheDocument()
  expect(screen.getByText('50% das 4 linhas')).toBeInTheDocument()
  expect(screen.getByText('Flag de qualidade, independente do estado do cálculo')).toBeInTheDocument()
  const legacy = screen.getByRole('region', { name: 'Dataset PM7' })
  expect(within(legacy).getByRole('row', { name: /Pending 1/ })).toBeInTheDocument()
  const store = screen.getByRole('region', { name: 'Semi-Imperium' })
  expect(within(store).getByRole('row', { name: /failed 2/ })).toBeInTheDocument()
})

test('informa quando a API está indisponível e permite reconectar', async () => {
  const fetchMock = mockApi({ '/api/health': [
    json({ detail: 'unavailable' }, 503),
    json({ status: 'ok', read_only: true }),
  ] })
  render(<App />)
  expect(await screen.findByText('API indisponível')).toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: 'Reconectar' }))
  expect(await screen.findByText('API ativa · somente leitura')).toBeInTheDocument()
  expect(fetchMock.mock.calls.filter(([path]) => path === '/api/health')).toHaveLength(2)
})

test('expõe falha da API e permite tentar novamente', async () => {
  let progressCalls = 0
  const fetchMock = vi.fn((input: RequestInfo | URL) => {
    if (String(input) === '/api/health') return Promise.resolve(json({ status: 'ok', read_only: true }))
    progressCalls += 1
    return progressCalls === 1 ? Promise.reject(new Error('offline')) : Promise.resolve(json(progress))
  })
  vi.stubGlobal('fetch', fetchMock)
  render(<App />)
  const alert = await screen.findByRole('alert')
  expect(alert).toHaveTextContent('Não foi possível conectar à API')
  fireEvent.click(within(alert).getByRole('button', { name: 'Tentar novamente' }))
  expect(await screen.findByText('Protocolo legado: crest_pm7')).toBeInTheDocument()
  expect(fetchMock.mock.calls.filter(([path]) => path === '/api/dataset/progress')).toHaveLength(2)
})

test('mostra arquivo ausente e estados vazios', async () => {
  mockApi({ '/api/dataset/progress': json({ detail: 'file not found', path: '/data/pm7.csv' }, 503) })
  const { unmount } = render(<App />)
  expect(await screen.findByRole('alert')).toHaveTextContent('/data/pm7.csv')
  unmount()
  mockApi({ '/api/dataset/progress': json({
    legacy: { ...progress.legacy, total: 0, status: {}, cbs_quality_flag: {} },
    semi_imperium: { ...progress.semi_imperium, exists: false, calculations: 0, states: {} },
  }) })
  render(<App />)
  expect(await screen.findByText('Nenhuma linha no dataset.')).toBeInTheDocument()
  expect(screen.getByText('Store ainda não existe.')).toBeInTheDocument()
})

test('permite filtrar validação, conflito e buscar ID', async () => {
  const fetchMock = mockApi()
  render(<App />)
  fireEvent.click(screen.getByRole('button', { name: 'Experimental' }))
  expect(await screen.findByText('Paridade CBS × experimental')).toBeInTheDocument()
  expect(screen.getByText('Erro por tamanho da molécula')).toBeInTheDocument()
  expect(screen.getByRole('img', { name: /Moléculas por família de fonte/ })).toHaveAccessibleName(/atct: 2; nist: 2/)
  expect(screen.getByRole('row', { name: /exp_2 cbs_2/ })).toHaveTextContent('só validação, conflito')
  fireEvent.click(screen.getByLabelText('Só validação'))
  expect(within(screen.getByLabelText('Estatísticas dos pares filtrados')).getByText('1')).toBeInTheDocument()
  expect(screen.queryByRole('row', { name: /exp_1 cbs_1/ })).not.toBeInTheDocument()
  fireEvent.change(screen.getByLabelText('Buscar ID, SMILES ou fonte'), { target: { value: 'cbs_1' } })
  expect(screen.getByText('Nenhum par corresponde aos filtros.')).toBeInTheDocument()
  expect(fetchMock).toHaveBeenCalledWith('/api/experimental/summary', expect.any(Object))
  expect(fetchMock).toHaveBeenCalledWith('/api/experimental/comparison', expect.any(Object))
})

test('mantém contagens se comparação falha e mostra resposta vazia', async () => {
  mockApi({ '/api/experimental/comparison': json({ detail: 'file not found', path: '/data/vs_cbs.csv' }, 503) })
  const { unmount } = render(<App />)
  fireEvent.click(screen.getByRole('button', { name: 'Experimental' }))
  expect(await screen.findByText('4', { selector: '.metric' })).toBeInTheDocument()
  expect(await screen.findByRole('alert')).toHaveTextContent('/data/vs_cbs.csv')
  unmount()
  mockApi({
    '/api/experimental/summary': json({ ...summary, rows: 0 }),
    '/api/experimental/comparison': json({ ...comparison, points: [] }),
  })
  render(<App />)
  fireEvent.click(screen.getByRole('button', { name: 'Experimental' }))
  expect(await screen.findByText('Nenhum valor experimental disponível.')).toBeInTheDocument()
  expect(await screen.findByText('Nenhum par CBS × experimental disponível.')).toBeInTheDocument()
})

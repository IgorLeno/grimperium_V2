import { cleanup, fireEvent, render, screen, within } from '@testing-library/react'
import type { ReactNode } from 'react'
import { afterEach, beforeEach, expect, test, vi } from 'vitest'
import { App } from './App'

vi.mock('recharts', () => {
  const container = ({ children }: { children?: ReactNode }) => <div>{children}</div>
  const empty = () => null
  return {
    Area: empty, AreaChart: container, Bar: empty, BarChart: container, CartesianGrid: empty, Cell: empty,
    ComposedChart: container, LabelList: empty, Legend: empty, Line: empty, Pie: container, PieChart: container,
    ReferenceLine: empty, ResponsiveContainer: container, Scatter: empty, ScatterChart: container,
    Tooltip: empty, XAxis: empty, YAxis: empty,
  }
})

const progress = {
  legacy: { protocol: 'crest_pm7', path: '/data/pm7.csv', total: 4,
    status: { OK: 2, Pending: 1, Failed: 1 }, cbs_quality_flag: { OK: 2, SUSPECT: 1, unset: 1 } },
  semi_imperium: { path: '/data/store', exists: true, molecules: 1, runs: 2,
    calculations: 3, states: { pending: 1, running: 0, verified: 0, unverified: 0, saddle: 0, failed: 2 } },
}
const summary = { path: '/data/exp.csv', rows: 4,
  by_source: { atct: 2, 'nist:gas': 2 }, by_source_family: { atct: 2, nist: 2 },
  flags: { conflict_flag: 1, validation_only: 1, low_priority_only: 0, stewart_reference: 0 },
  uncertainty_kind: { reported: 3, none: 1 } }
const comparison = { path: '/data/vs_cbs.csv', units: 'kcal/mol; cbs_minus_exp = H298_cbs - H298_exp',
  stats: { n: 2, mae: 2, bias: 2, rmse: 2.24 }, points: [
    { exp_id: 'exp_1', mol_id: 'cbs_1', smiles: 'C', nheavy: 1, source: 'atct', H298_exp: -17, H298_cbs: -16, cbs_minus_exp: 1, uncertainty: .1, conflict_flag: false, low_priority_only: false, validation_only: false },
    { exp_id: 'exp_2', mol_id: 'cbs_2', smiles: 'CCO', nheavy: 3, source: 'nist:gas', H298_exp: -56, H298_cbs: -53, cbs_minus_exp: 3, uncertainty: null, conflict_flag: true, low_priority_only: false, validation_only: true },
  ] }

function histogram(column: string, label: string, n: number) {
  return { column, label, unit: 'kcal/mol', n, below: 0, above: 0, mean: -5, median: -5, p05: -10, p95: 0,
    minimum: -10, maximum: 0, bins: [{ start: -10, end: -5, count: Math.ceil(n / 2) }, { start: -5, end: 0, count: Math.floor(n / 2) }] }
}

function profile(key: string, rows: number, extra: Record<string, unknown> = {}) {
  return { key, title: key, path: `/data/${key}.csv`, rows, unique_smiles: rows, nheavy: { 1: 1, 3: rows - 1 },
    composition: { CH: 1, CHO: rows - 1, CHN: 0, CHNO: 0 }, categories: {}, histograms: [], by_nheavy: [],
    pairs: [], pairs_total: 0, missing_columns: [], manifest: null, ...extra }
}

const expProfile = profile('exp', 4, { histograms: [histogram('H298_exp', 'Entalpia', 4)] })
const cbsProfile = profile('cbs', 10, {
  categories: { multiplicity: { 1: 9, 2: 1 }, charge: { 0: 10 }, smiles_count: { 1: 8, 2: 2 } },
  histograms: [histogram('H298_cbs', 'Entalpia CBS', 10), histogram('cbs_b3', 'Correção', 10)],
  manifest: { source_citation: 'Dobbelaere et al. 2024' },
})
const pm7Profile = profile('pm7', 4, {
  categories: { status: { OK: 2, Pending: 1, Failed: 1 }, cbs_quality_flag: { OK: 2, unset: 2 } },
  histograms: [histogram('delta_cbs_pm7', 'Erro', 2)],
  by_nheavy: [{ nheavy: 1, n: 1, mean: -2, mean_abs: 2 }],
  pairs: [{ id: 'a', x: -17, y: -15, nheavy: 1 }], pairs_total: 1,
})
const rows = (key: string) => ({ key, columns: ['mol_id', 'smiles', 'H298_cbs'], total: 2, matched: 2, offset: 0, limit: 25,
  rows: [{ mol_id: `${key}_row_1`, smiles: 'C', H298_cbs: -17.8 }, { mol_id: `${key}_row_2`, smiles: 'CCO', H298_cbs: -56.2 }] })

const catalog = {
  dispatch_enabled: false,
  properties: [{ key: 'standard_enthalpy_of_formation', label: 'Entalpia de formação (ΔHf° 298 K)', unit: 'kcal/mol', description: 'ΔHf', available: true }],
  input_modes: [
    { key: 'smiles', label: 'SMILES', role: 'input', description: 'Uma por linha', executable: null, available: null, experimental: false, parameters: [] },
    { key: 'name', label: 'Nome da molécula (PubChem)', role: 'input', description: 'PubChem', executable: null, available: null, experimental: false, parameters: [] },
  ],
  programs: ['rdkit', 'crest', 'confpass', 'mopac'].map((key) => ({
    key, label: key.toUpperCase(), role: 'step', description: key, executable: key === 'crest' || key === 'mopac' ? key : null,
    available: key === 'crest' ? false : null, experimental: key === 'confpass',
    parameters: key === 'mopac' ? [{ key: 'hamiltonian', label: 'Hamiltoniano', kind: 'choice', default: 'PM7', choices: ['AM1', 'PM3', 'PM7'], unit: '', help: '', minimum: null, maximum: null }] : [],
  })),
  outputs: [{ key: 'spreadsheet', label: 'Planilha', role: 'output', description: 'CSV', executable: null, available: null, experimental: false, parameters: [] }],
}

function json(data: unknown, status = 200): Response {
  return new Response(JSON.stringify(data), { status, headers: { 'Content-Type': 'application/json' } })
}

function mockApi(overrides: Record<string, () => Response> = {}) {
  const routes: Record<string, () => Response> = {
    '/api/health': () => json({ status: 'ok', read_only: true }),
    '/api/dataset/progress': () => json(progress),
    '/api/experimental/summary': () => json(summary),
    '/api/experimental/comparison': () => json(comparison),
    '/api/tables/exp/profile': () => json(expProfile),
    '/api/tables/cbs/profile': () => json(cbsProfile),
    '/api/tables/pm7/profile': () => json(pm7Profile),
    '/api/routines/catalog': () => json(catalog),
    ...overrides,
  }
  const fetchMock = vi.fn((input: RequestInfo | URL) => {
    const path = String(input)
    const table = path.match(/^\/api\/tables\/(\w+)\/rows/)
    if (table) return Promise.resolve(json(rows(table[1])))
    const route = routes[path]
    return route ? Promise.resolve(route()) : Promise.resolve(json({ detail: 'not mocked' }, 404))
  })
  vi.stubGlobal('fetch', fetchMock)
  return fetchMock
}

function openTab(name: RegExp) {
  const menu = within(screen.getByRole('navigation', { name: 'Telas' }))
  const top = menu.queryByRole('button', { name })
  if (top) { fireEvent.click(top); return }
  if (!screen.queryByRole('navigation', { name: 'Submenu de Dados' })) fireEvent.click(menu.getByRole('button', { name: /Dados/ }))
  fireEvent.click(within(screen.getByRole('navigation', { name: 'Submenu de Dados' })).getByRole('button', { name }))
}

/** O app abre no hub Delta Chemistry; as telas de dados vivem no módulo Grimperium. */
function renderGrimperium() {
  const view = render(<App />)
  fireEvent.click(within(screen.getByRole('navigation', { name: 'Módulos de propriedade' })).getByRole('button', { name: /Grimperium/ }))
  return view
}

beforeEach(() => { window.localStorage.clear() })
afterEach(() => { cleanup(); vi.unstubAllGlobals() })

test('resumo escreve o estado do projeto com números reais', async () => {
  mockApi()
  renderGrimperium()
  expect(await screen.findByText('Protocolo legado: crest_pm7')).toBeInTheDocument()
  expect(screen.getByText('API ativa · somente leitura')).toBeInTheDocument()
  const text = screen.getByRole('region', { name: 'Resumo em texto' })
  await within(text).findByText('10')
  expect(text).toHaveTextContent(/reúne 10 geometrias/)
  expect(text).toHaveTextContent(/tem 4 moléculas medidas/)
  expect(text).toHaveTextContent(/das quais 2 também têm valor CBS/)
  expect(text).toHaveTextContent(/concluiu 2 de 4 moléculas \(50%\)/)
  expect(screen.getByText('Flag de qualidade, independente do estado do cálculo')).toBeInTheDocument()
  const legacy = screen.getByRole('region', { name: 'Dataset PM7' })
  expect(within(legacy).getByRole('img', { name: /Estados do dataset PM7/ })).toHaveAccessibleName(/Na fila \(Pending\): 1/)
  const store = screen.getByRole('region', { name: 'Semi-Imperium' })
  expect(within(store).getByRole('img', { name: /Estado mais recente/ })).toHaveAccessibleName(/Falhou: 2/)
})

test('informa quando a API está indisponível e permite reconectar', async () => {
  let calls = 0
  const fetchMock = mockApi({ '/api/health': () => (calls++ === 0 ? json({ detail: 'unavailable' }, 503) : json({ status: 'ok', read_only: true })) })
  render(<App />)
  expect(await screen.findByText('API indisponível')).toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: 'Reconectar' }))
  expect(await screen.findByText('API ativa · somente leitura')).toBeInTheDocument()
  expect(fetchMock.mock.calls.filter(([path]) => path === '/api/health')).toHaveLength(2)
})

test('mostra estados vazios do dataset legado e do store', async () => {
  mockApi({ '/api/dataset/progress': () => json({
    legacy: { ...progress.legacy, total: 0, status: {}, cbs_quality_flag: {} },
    semi_imperium: { ...progress.semi_imperium, exists: false, calculations: 0, states: {} },
  }) })
  renderGrimperium()
  expect(await screen.findByText('Nenhuma linha no dataset.')).toBeInTheDocument()
  expect(screen.getByText(/Store ainda não existe/)).toBeInTheDocument()
})

test('aba experimental filtra validação, conflito e busca ID', async () => {
  const fetchMock = mockApi()
  renderGrimperium()
  openTab(/Experimental/)
  expect(await screen.findByText('Paridade CBS × experimental')).toBeInTheDocument()
  expect(screen.getByText('Erro por tamanho da molécula')).toBeInTheDocument()
  expect(screen.getByText('O que é a tabela experimental?')).toBeInTheDocument()
  expect(screen.getByRole('img', { name: /Moléculas por família de fonte/ })).toHaveAccessibleName(/atct: 2; nist: 2/)
  expect(screen.getByRole('row', { name: /exp_2 cbs_2/ })).toHaveTextContent('só validação, conflito')
  fireEvent.click(screen.getByLabelText('Só validação'))
  const stats = screen.getByLabelText('Estatísticas dos pares filtrados')
  expect(within(stats).getByText('Pares no recorte').parentElement?.parentElement).toHaveTextContent('1')
  expect(screen.queryByRole('row', { name: /exp_1 cbs_1/ })).not.toBeInTheDocument()
  fireEvent.change(screen.getByLabelText('Buscar ID, SMILES ou fonte'), { target: { value: 'cbs_1' } })
  expect(screen.getByText('Nenhum par corresponde aos filtros.')).toBeInTheDocument()
  expect(await screen.findByRole('region', { name: 'Valores experimentais' })).toBeInTheDocument()
  expect(fetchMock).toHaveBeenCalledWith('/api/experimental/comparison', expect.any(Object))
})

test('aba experimental mantém contagens se a comparação falha', async () => {
  mockApi({ '/api/experimental/comparison': () => json({ detail: 'file not found', path: '/data/vs_cbs.csv' }, 503) })
  renderGrimperium()
  openTab(/Experimental/)
  expect(await screen.findByText('4', { selector: '.metric' })).toBeInTheDocument()
  const alerts = await screen.findAllByRole('alert')
  expect(alerts.some((alert) => alert.textContent?.includes('/data/vs_cbs.csv'))).toBe(true)
})

test('aba CBS explica a tabela e permite ordenar a tabela completa', async () => {
  const fetchMock = mockApi()
  renderGrimperium()
  openTab(/CBS/)
  expect(await screen.findByText('O que é a tabela CBS?')).toBeInTheDocument()
  expect(await screen.findByText('Dobbelaere et al. 2024')).toBeInTheDocument()
  expect(await screen.findByText('cbs_row_1')).toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: /ΔHf CBS/ }))
  await screen.findByText('cbs_row_1')
  expect(fetchMock.mock.calls.some(([path]) => String(path).includes('/api/tables/cbs/rows') && String(path).includes('sort=H298_cbs'))).toBe(true)
})

test('aba CBS explica quando o CSV não está no computador', async () => {
  mockApi({ '/api/tables/cbs/profile': () => json({ detail: 'file not found', path: '/data/thermo_cbs_chon.csv' }, 503) })
  renderGrimperium()
  openTab(/CBS/)
  expect(await screen.findByRole('region', { name: 'Arquivo não encontrado' })).toHaveTextContent('data/thermo_cbs_chon.csv')
})

test('aba PM7 mostra andamento e erro do cálculo rápido', async () => {
  mockApi()
  renderGrimperium()
  openTab(/PM7/)
  expect(await screen.findByText('O que é a tabela PM7?')).toBeInTheDocument()
  expect(await screen.findByRole('img', { name: /Moléculas concluídas: 50%/ })).toBeInTheDocument()
  expect(screen.getByText('Paridade CBS × PM7')).toBeInTheDocument()
  expect(screen.getByText('Protocolo legado: crest_pm7')).toBeInTheDocument()
})

test('rotinas: monta, valida e salva uma rotina sem disparar', async () => {
  mockApi()
  renderGrimperium()
  openTab(/Calcular/)
  fireEvent.click((await screen.findAllByRole('button', { name: /Criar nova rotina/ }))[0])
  fireEvent.change(screen.getByLabelText('Nome da rotina'), { target: { value: 'Teste álcoois' } })
  fireEvent.click(screen.getByRole('button', { name: /Pipeline/ }))
  fireEvent.change(screen.getByLabelText(/Estruturas/), { target: { value: 'CCO\nCCCO' } })
  expect(screen.getByText('2 molécula(s) informada(s).')).toBeInTheDocument()
  fireEvent.click(screen.getByRole('checkbox', { name: 'MOPAC' }))
  expect(screen.getByText('O MOPAC é obrigatório para calcular a entalpia de formação.')).toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: /Direto no PM7/ }))
  expect(screen.getByText('✓ Ordem válida para o pipeline Semi-Imperium.')).toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: /Revisão/ }))
  expect(screen.getByRole('button', { name: /Disparar/ })).toBeDisabled()
  fireEvent.click(screen.getByRole('button', { name: 'Salvar rotina' }))
  expect(await screen.findByText('Rotina “Teste álcoois” salva.')).toBeInTheDocument()
  const menu = screen.getByRole('complementary', { name: 'Menu de rotinas' })
  expect(within(menu).getByText('Teste álcoois')).toBeInTheDocument()
  expect(JSON.parse(window.localStorage.getItem('grimperium.routines.v1') ?? '[]')[0].steps).toEqual(['rdkit', 'mopac'])
})

test('abre no hub Delta Chemistry: Rotinas, Resumos e Sobre, sem as abas do Grimperium', async () => {
  mockApi()
  render(<App />)
  const menu = within(screen.getByRole('navigation', { name: 'Telas' }))
  expect(menu.getAllByRole('button').map((button) => button.textContent)).toEqual(['Rotinas', 'Resumos', 'Sobre'])
  expect(menu.getByRole('button', { name: 'Rotinas' })).toHaveAttribute('aria-current', 'page')
  expect(await screen.findByText('Rotinas · Painel de controle')).toBeInTheDocument()
  expect(menu.queryByRole('button', { name: /Dados|Calcular|Treinar/ })).not.toBeInTheDocument()
  fireEvent.click(menu.getByRole('button', { name: 'Resumos' }))
  expect(await screen.findByText('Nenhum resultado ainda')).toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: 'Ir para Rotinas' }))
  expect(await screen.findByText('Rotinas · Painel de controle')).toBeInTheDocument()
  fireEvent.click(menu.getByRole('button', { name: 'Sobre' }))
  expect(await screen.findByText('Rapidez com precisão, em três passos')).toBeInTheDocument()
  expect(screen.getByText(/Hoje só existe o módulo de entalpia de formação/)).toBeInTheDocument()
})

test('barra lateral: hub Delta Chemistry e módulo Grimperium, que traz as telas de dados', async () => {
  mockApi()
  render(<App />)
  const sidebar = screen.getByRole('complementary', { name: 'Navegação lateral' })
  const hub = within(within(sidebar).getByRole('navigation', { name: 'Hub' })).getByRole('button', { name: /Delta Chemistry/ })
  expect(hub).toHaveAttribute('aria-current', 'true')
  const module = within(within(sidebar).getByRole('navigation', { name: 'Módulos de propriedade' })).getByRole('button', { name: /Grimperium/ })
  expect(module).toHaveTextContent('Entalpia de formação (estado padrão)')
  expect(module).not.toHaveAttribute('aria-current')
  expect(within(sidebar).getAllByRole('button')).toHaveLength(2)
  fireEvent.click(module)
  expect(module).toHaveAttribute('aria-current', 'true')
  expect(hub).not.toHaveAttribute('aria-current')
  expect(screen.getByLabelText('Local atual')).toHaveTextContent('Delta Chemistry/Grimperium')
  openTab(/CBS/)
  expect(await screen.findByText('O que é a tabela CBS?')).toBeInTheDocument()
  fireEvent.click(module)
  expect(within(screen.getByRole('navigation', { name: 'Submenu de Dados' })).getByRole('button', { name: 'CBS' })).toHaveAttribute('aria-current', 'page')
  fireEvent.click(hub)
  expect(within(screen.getByRole('navigation', { name: 'Telas' })).getAllByRole('button').map((button) => button.textContent)).toEqual(['Rotinas', 'Resumos', 'Sobre'])
})

test('módulo Grimperium: Resumo, Dados com submenu, Calcular e Treinar', async () => {
  mockApi()
  renderGrimperium()
  const menu = within(screen.getByRole('navigation', { name: 'Telas' }))
  expect(menu.getAllByRole('button').map((button) => button.textContent?.replace('▾', ''))).toEqual(['Resumo', 'Dados', 'Calcular', 'Treinar'])
  expect(screen.queryByRole('navigation', { name: 'Submenu de Dados' })).not.toBeInTheDocument()
  fireEvent.click(menu.getByRole('button', { name: /Dados/ }))
  const sub = within(screen.getByRole('navigation', { name: 'Submenu de Dados' }))
  expect(sub.getAllByRole('button').map((button) => button.textContent)).toEqual(['Experimental', 'CBS', 'PM7'])
  expect(await screen.findByText('Paridade CBS × experimental')).toBeInTheDocument()
  fireEvent.click(sub.getByRole('button', { name: 'PM7' }))
  expect(await screen.findByText('O que é a tabela PM7?')).toBeInTheDocument()
  expect(menu.getByRole('button', { name: /Dados/ })).toHaveAttribute('aria-current', 'page')
  fireEvent.click(menu.getByRole('button', { name: 'Treinar' }))
  expect(await screen.findByText('Cinco etapas, sem vazamento de dados')).toBeInTheDocument()
  expect(screen.queryByRole('navigation', { name: 'Submenu de Dados' })).not.toBeInTheDocument()
})

test('cabeçalho só tem a marca; o botão do menu lateral fica abaixo dele', () => {
  mockApi()
  render(<App />)
  const header = screen.getByRole('banner')
  expect(within(header).getByText('Delta Chemistry')).toBeInTheDocument()
  const toggle = screen.getByRole('button', { name: 'Ocultar menu lateral' })
  expect(header).not.toContainElement(toggle)
  expect(header.compareDocumentPosition(toggle) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy()
})

test('menu lateral pode ser ocultado e a escolha é lembrada', () => {
  mockApi()
  const { unmount } = render(<App />)
  expect(screen.getByRole('button', { name: 'Ocultar menu lateral' })).toHaveAttribute('aria-expanded', 'true')
  fireEvent.click(screen.getByRole('button', { name: 'Ocultar menu lateral' }))
  expect(screen.queryByRole('complementary', { name: 'Navegação lateral' })).not.toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Mostrar menu lateral' })).toHaveAttribute('aria-expanded', 'false')
  expect(window.localStorage.getItem('grimperium.sidebar.open')).toBe('false')
  unmount()
  render(<App />)
  expect(screen.queryByRole('complementary', { name: 'Navegação lateral' })).not.toBeInTheDocument()
  fireEvent.click(screen.getByRole('button', { name: 'Mostrar menu lateral' }))
  expect(screen.getByRole('complementary', { name: 'Navegação lateral' })).toBeInTheDocument()
})

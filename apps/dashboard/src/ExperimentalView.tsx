import { useMemo, useState } from 'react'
import {
  CartesianGrid,
  Legend,
  ReferenceLine,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis,
} from 'recharts'
import type { ComparisonPoint, ExperimentalComparison, ExperimentalSummary } from './api'
import { CountTable, ErrorNotice, Section } from './Resource'
import { useApi } from './useApi'

const PAGE_SIZE = 25

type Category = 'normal' | 'validation' | 'conflict' | 'both'

const categories: { key: Category; label: string; color: string }[] = [
  { key: 'normal', label: 'Demais pares', color: '#3977aa' },
  { key: 'validation', label: 'Só validação', color: '#ad6f1f' },
  { key: 'conflict', label: 'Conflito', color: '#b34048' },
  { key: 'both', label: 'Validação + conflito', color: '#7543a2' },
]

function category(point: ComparisonPoint): Category {
  if (point.validation_only && point.conflict_flag) return 'both'
  if (point.validation_only) return 'validation'
  if (point.conflict_flag) return 'conflict'
  return 'normal'
}

function format(value: number | null): string {
  return value === null ? '—' : value.toLocaleString('pt-BR', { maximumFractionDigits: 2 })
}

function FilteredStats({ points }: { points: ComparisonPoint[] }) {
  const count = points.length
  const bias = count ? points.reduce((sum, point) => sum + point.cbs_minus_exp, 0) / count : null
  const mae = count ? points.reduce((sum, point) => sum + Math.abs(point.cbs_minus_exp), 0) / count : null
  const rmse = count ? Math.sqrt(points.reduce((sum, point) => sum + point.cbs_minus_exp ** 2, 0) / count) : null
  return (
    <div className="metrics" aria-label="Estatísticas dos pares filtrados">
      <p><strong>{count.toLocaleString('pt-BR')}</strong> pares filtrados</p>
      <p><strong>{format(bias)}</strong> viés CBS − exp</p>
      <p><strong>{format(mae)}</strong> MAE</p>
      <p><strong>{format(rmse)}</strong> RMSE</p>
    </div>
  )
}

function PairCharts({ points }: { points: ComparisonPoint[] }) {
  const values = points.flatMap((point) => [point.H298_exp, point.H298_cbs])
  const lower = Math.floor(Math.min(...values))
  const upper = Math.ceil(Math.max(...values))
  const parityDomain = lower === upper ? [lower - 1, upper + 1] : [lower, upper]
  const groups = categories.map((item) => ({
    ...item,
    data: points.filter((point) => category(point) === item.key),
  }))
  return (
    <div className="columns charts">
      <Section title="Paridade CBS × experimental">
        <p>ΔHf em kcal/mol; linha tracejada: igualdade.</p>
        <div className="chart" role="img" aria-label="Dispersão de entalpia CBS versus experimental">
          <ResponsiveContainer width="100%" height={310}>
            <ScatterChart margin={{ top: 12, right: 18, bottom: 18, left: 14 }}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis type="number" dataKey="H298_exp" name="Experimental" domain={parityDomain} label={{ value: 'Experimental', position: 'insideBottom', offset: -10 }} />
              <YAxis type="number" dataKey="H298_cbs" name="CBS" domain={parityDomain} label={{ value: 'CBS', angle: -90, position: 'insideLeft' }} />
              <Tooltip cursor={{ strokeDasharray: '3 3' }} />
              <Legend verticalAlign="top" />
              <ReferenceLine segment={[{ x: parityDomain[0], y: parityDomain[0] }, { x: parityDomain[1], y: parityDomain[1] }]} stroke="#777" strokeDasharray="4 4" />
              {groups.filter((group) => group.data.length).map((group) => (
                <Scatter key={group.key} name={group.label} data={group.data} fill={group.color} />
              ))}
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      </Section>
      <Section title="Erro por tamanho da molécula">
        <p>CBS − experimental (kcal/mol) × número de átomos pesados.</p>
        <div className="chart" role="img" aria-label="Dispersão de erro CBS menos experimental por número de átomos pesados">
          <ResponsiveContainer width="100%" height={310}>
            <ScatterChart margin={{ top: 12, right: 18, bottom: 18, left: 14 }}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis type="number" dataKey="nheavy" name="Átomos pesados" allowDecimals={false} label={{ value: 'Átomos pesados', position: 'insideBottom', offset: -10 }} />
              <YAxis type="number" dataKey="cbs_minus_exp" name="CBS − exp" label={{ value: 'CBS − exp', angle: -90, position: 'insideLeft' }} />
              <Tooltip cursor={{ strokeDasharray: '3 3' }} />
              <Legend verticalAlign="top" />
              <ReferenceLine y={0} stroke="#777" strokeDasharray="4 4" />
              {groups.filter((group) => group.data.length).map((group) => (
                <Scatter key={group.key} name={group.label} data={group.data} fill={group.color} />
              ))}
            </ScatterChart>
          </ResponsiveContainer>
        </div>
      </Section>
    </div>
  )
}

function ComparisonExplorer({ comparison }: { comparison: ExperimentalComparison }) {
  const [query, setQuery] = useState('')
  const [source, setSource] = useState('all')
  const [validationOnly, setValidationOnly] = useState(false)
  const [conflictOnly, setConflictOnly] = useState(false)
  const [page, setPage] = useState(0)
  const sources = useMemo(
    () => [...new Set(comparison.points.map((point) => point.source))].sort(),
    [comparison.points],
  )
  const filtered = useMemo(() => {
    const term = query.trim().toLowerCase()
    return comparison.points
      .filter((point) =>
        (source === 'all' || point.source === source) &&
        (!validationOnly || point.validation_only) &&
        (!conflictOnly || point.conflict_flag) &&
        (!term || [point.exp_id, point.mol_id, point.smiles, point.source].some((value) => value.toLowerCase().includes(term))),
      )
      .sort((a, b) => Math.abs(b.cbs_minus_exp) - Math.abs(a.cbs_minus_exp))
  }, [comparison.points, query, source, validationOnly, conflictOnly])
  const maxPage = Math.max(0, Math.ceil(filtered.length / PAGE_SIZE) - 1)
  const visible = filtered.slice(Math.min(page, maxPage) * PAGE_SIZE, (Math.min(page, maxPage) + 1) * PAGE_SIZE)
  const resetPage = () => setPage(0)

  return (
    <>
      <Section title="Explorar pares CBS × experimental">
        <p>Valores NIST marcados como “só validação” não são alvos de treino. Pares ordenados por maior |CBS − exp|.</p>
        <div className="filters">
          <label>Buscar ID, SMILES ou fonte
            <input value={query} onChange={(event) => { setQuery(event.target.value); resetPage() }} placeholder="Ex.: cbs_01102" />
          </label>
          <label>Fonte
            <select value={source} onChange={(event) => { setSource(event.target.value); resetPage() }}>
              <option value="all">Todas</option>
              {sources.map((item) => <option key={item} value={item}>{item}</option>)}
            </select>
          </label>
          <label className="checkbox"><input type="checkbox" checked={validationOnly} onChange={(event) => { setValidationOnly(event.target.checked); resetPage() }} /> Só validação</label>
          <label className="checkbox"><input type="checkbox" checked={conflictOnly} onChange={(event) => { setConflictOnly(event.target.checked); resetPage() }} /> Só conflitos</label>
        </div>
        <FilteredStats points={filtered} />
        {filtered.length === 0 ? (
          <p>Nenhum par corresponde aos filtros.</p>
        ) : (
          <>
            <PairCharts points={filtered} />
            <div className="table-wrap">
              <table>
                <thead><tr><th scope="col">ID exp.</th><th scope="col">ID CBS</th><th scope="col">SMILES</th><th scope="col">Fonte</th><th scope="col">Pesados</th><th scope="col">Exp.</th><th scope="col">CBS</th><th scope="col">CBS − exp</th><th scope="col">Flags</th></tr></thead>
                <tbody>{visible.map((point) => (
                  <tr key={point.exp_id}>
                    <td>{point.exp_id}</td><td>{point.mol_id}</td><td className="smiles">{point.smiles}</td><td>{point.source}</td>
                    <td>{point.nheavy}</td><td>{format(point.H298_exp)}</td><td>{format(point.H298_cbs)}</td><td>{format(point.cbs_minus_exp)}</td>
                    <td>{[point.validation_only && 'só validação', point.conflict_flag && 'conflito', point.low_priority_only && 'baixa prioridade'].filter(Boolean).join(', ') || '—'}</td>
                  </tr>
                ))}</tbody>
              </table>
            </div>
            <div className="pagination">
              <span>Mostrando {Math.min(page, maxPage) * PAGE_SIZE + 1}–{Math.min((Math.min(page, maxPage) + 1) * PAGE_SIZE, filtered.length)} de {filtered.length.toLocaleString('pt-BR')}</span>
              <button type="button" disabled={page === 0} onClick={() => setPage((value) => value - 1)}>Anterior</button>
              <button type="button" disabled={page >= maxPage} onClick={() => setPage((value) => value + 1)}>Próxima</button>
            </div>
          </>
        )}
      </Section>
    </>
  )
}

export function ExperimentalView() {
  const summary = useApi<ExperimentalSummary>('/api/experimental/summary')
  const comparison = useApi<ExperimentalComparison>('/api/experimental/comparison')
  return (
    <div className="view">
      <div className="view-heading">
        <div><h1>Experimental</h1><p>Fontes de ΔHf gasoso e diferenças contra a referência CBS.</p></div>
        <button type="button" onClick={() => { summary.reload(); comparison.reload() }}>Atualizar dados</button>
      </div>
      <Section title="Cobertura e proveniência">
        {summary.state.status === 'loading' && <p role="status">Carregando fontes…</p>}
        {summary.state.status === 'error' && <ErrorNotice message={summary.state.message} retry={summary.reload} />}
        {summary.state.status === 'ready' && (summary.state.data.rows === 0 ? <p>Nenhum valor experimental disponível.</p> : (
          <>
            <p className="metric">{summary.state.data.rows.toLocaleString('pt-BR')} <small>moléculas</small></p>
            <div className="columns">
              <div><h3>Por família</h3><CountTable counts={summary.state.data.by_source_family} label="Família" /></div>
              <div><h3>Por fonte</h3><CountTable counts={summary.state.data.by_source} label="Fonte" /></div>
              <div><h3>Flags</h3><CountTable counts={summary.state.data.flags} label="Flag" /><h3>Incerteza</h3><CountTable counts={summary.state.data.uncertainty_kind} label="Tipo" /></div>
            </div>
          </>
        ))}
      </Section>
      {comparison.state.status === 'loading' && <p role="status">Carregando comparação CBS…</p>}
      {comparison.state.status === 'error' && <Section title="Comparação CBS"><ErrorNotice message={comparison.state.message} retry={comparison.reload} /></Section>}
      {comparison.state.status === 'ready' && (comparison.state.data.points.length === 0 ? (
        <Section title="Comparação CBS"><p>Nenhum par CBS × experimental disponível.</p></Section>
      ) : <ComparisonExplorer comparison={comparison.state.data} />)}
    </div>
  )
}

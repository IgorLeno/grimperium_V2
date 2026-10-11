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
import { MetricCard, MoleculeIllustration, RankedBars, StoryNote } from './StoryVisuals'
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
    <div className="metric-grid metric-grid--four" aria-label="Estatísticas dos pares filtrados">
      <MetricCard label="Pares no recorte" value={count.toLocaleString('pt-BR')} note="Moléculas com CBS e medida experimental" accent />
      <MetricCard label="Diferença média" value={format(bias)} note="CBS − experimento · kcal/mol" />
      <MetricCard label="Erro absoluto médio" value={format(mae)} note="Distância média, sem sinal · kcal/mol" />
      <MetricCard label="RMSE" value={format(rmse)} note="Raiz do erro quadrático médio; destaca diferenças grandes · kcal/mol" />
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
        <div className="section-heading"><div><span className="eyebrow eyebrow--dark">Etapa 02 / comparação</span><h3>Quando cálculo e medida se encontram</h3><p>Um par é uma molécula que tem ambos os valores. A distância entre eles é medida em {comparison.units}.</p></div></div>
        <div className="editorial-callout"><span>Como ler</span><p><strong>Quanto mais perto de zero, menor a diferença.</strong> O sinal indica a direção: positivo significa CBS acima do experimental. Não é uma medida de causa.</p></div>
        <div className="filter-heading"><div><span className="eyebrow eyebrow--dark">Explore o conjunto</span><h3>Escolha o recorte</h3></div><p>Indicadores, gráficos e tabela respondem aos mesmos filtros.</p></div>
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
            <div className="table-heading"><div><span className="eyebrow eyebrow--dark">Os dados em detalhe</span><h3>Molécula por molécula</h3><p>Ordenado pela maior diferença absoluta entre CBS e experimento. Valores em {comparison.units}.</p></div></div>
            <div className="table-wrap">
              <table className="pair-table">
                <caption className="sr-only">Pares CBS e experimental filtrados, ordenados pela diferença absoluta</caption>
                <thead><tr><th scope="col">ID exp.</th><th scope="col">ID CBS</th><th scope="col">SMILES</th><th scope="col">Fonte</th><th scope="col">Pesados</th><th scope="col">Exp.</th><th scope="col">CBS</th><th scope="col">CBS − exp</th><th scope="col">Flags</th></tr></thead>
                <tbody>{visible.map((point) => (
                  <tr key={point.exp_id}>
                    <td className="id-cell">{point.exp_id}</td><td className="id-cell">{point.mol_id}</td><td className="smiles">{point.smiles}</td><td>{point.source}</td>
                    <td className="num-cell">{point.nheavy}</td><td className="num-cell">{format(point.H298_exp)}</td><td className="num-cell">{format(point.H298_cbs)}</td><td className="num-cell"><span className={Math.abs(point.cbs_minus_exp) >= 10 ? 'difference difference--high' : 'difference'}>{format(point.cbs_minus_exp)}</span></td>
                    <td className="flags-cell">{[point.validation_only && 'só validação', point.conflict_flag && 'conflito', point.low_priority_only && 'baixa prioridade'].filter(Boolean).join(', ') || '—'}</td>
                  </tr>
                ))}</tbody>
              </table>
            </div>
            <div className="pagination">
              <span>Mostrando {Math.min(page, maxPage) * PAGE_SIZE + 1}–{Math.min((Math.min(page, maxPage) + 1) * PAGE_SIZE, filtered.length)} de {filtered.length.toLocaleString('pt-BR')}</span>
              <button type="button" disabled={page === 0} onClick={() => setPage((value) => value - 1)}>Anterior</button>
              <button type="button" disabled={page >= maxPage} onClick={() => setPage((value) => value + 1)}>Próxima</button>
            </div>
            <p className="table-footnote">Valores “só validação” não são alvos de treino. Flags podem coexistir na mesma molécula.</p>
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
      <div className="hero hero--experimental">
        <div className="hero__copy">
          <span className="eyebrow">02 / Referência experimental</span>
          <h1>O cálculo encontra<br /><em>o mundo medido.</em></h1>
          <p>Valores experimentais ajudam a entender a diferença entre uma previsão computacional e uma medida de laboratório.</p>
          <button type="button" className="refresh-button" onClick={() => { summary.reload(); comparison.reload() }}>↻ <span>Atualizar dados</span></button>
        </div>
        <div className="hero__art"><MoleculeIllustration /><span>Uma representação visual da rede molecular</span></div>
      </div>
      <Section title="Cobertura e proveniência">
        {summary.state.status === 'loading' && <p role="status">Carregando fontes…</p>}
        {summary.state.status === 'error' && <ErrorNotice message={summary.state.message} retry={summary.reload} />}
        {summary.state.status === 'ready' && (summary.state.data.rows === 0 ? <p>Nenhum valor experimental disponível.</p> : (
          <>
            <div className="section-heading"><div><span className="eyebrow eyebrow--dark">Etapa 01 / cobertura</span><h3>De onde vêm as medidas?</h3><p>Cada molécula tem uma fonte selecionada. As barras comparam famílias de fontes no conjunto experimental; os percentuais usam todas as {summary.state.data.rows.toLocaleString('pt-BR')} moléculas como base.</p></div></div>
            <div className="coverage-layout">
              <div className="coverage-total"><span>Base experimental</span><p className="metric">{summary.state.data.rows.toLocaleString('pt-BR')}</p><strong>moléculas com valor experimental</strong><p>ΔHf° gasoso a 298,15 K</p></div>
              <div className="coverage-chart"><h3>Famílias de fontes</h3><p className="subtle">Participação sobre {summary.state.data.rows.toLocaleString('pt-BR')} moléculas</p><RankedBars counts={summary.state.data.by_source_family} total={summary.state.data.rows} label="Moléculas por família de fonte" /></div>
            </div>
            <div className="detail-grid detail-grid--three">
              <div><h3>Fonte selecionada</h3><p className="subtle">Uma fonte por molécula.</p><CountTable counts={summary.state.data.by_source} label="Fonte" /></div>
              <div><h3>Observações dos dados</h3><p className="subtle">Flags podem se sobrepor; não são partes de um total.</p><CountTable counts={summary.state.data.flags} label="Flag" /></div>
              <div><h3>Incerteza informada</h3><p className="subtle">Como a incerteza foi descrita na fonte.</p><CountTable counts={summary.state.data.uncertainty_kind} label="Tipo" /></div>
            </div>
          </>
        ))}
      </Section>
      <div className="story-grid">
        <StoryNote number="01" title="Medida e cálculo são coisas diferentes">O valor experimental vem de uma fonte de referência; CBS é um cálculo. A comparação mostra a distância entre os dois.</StoryNote>
        <StoryNote number="02" title="Nem todo valor serve para treino">Registros marcados “só validação” permanecem úteis para comparação, mas não são alvos de treinamento do modelo.</StoryNote>
      </div>
      {comparison.state.status === 'loading' && <p role="status">Carregando comparação CBS…</p>}
      {comparison.state.status === 'error' && <Section title="Comparação CBS"><ErrorNotice message={comparison.state.message} retry={comparison.reload} /></Section>}
      {comparison.state.status === 'ready' && (comparison.state.data.points.length === 0 ? (
        <Section title="Comparação CBS"><p>Nenhum par CBS × experimental disponível.</p></Section>
      ) : <ComparisonExplorer comparison={comparison.state.data} />)}
    </div>
  )
}

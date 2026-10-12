import { useMemo, useState } from 'react'
import type { ComparisonPoint, DatasetProfile, ExperimentalComparison, ExperimentalSummary } from './api'
import {
  CumulativeChart, DonutChart, GroupErrorChart, HistogramChart, LegendRow, ParityChart, SizeErrorChart,
} from './Charts'
import type { ParityDatum } from './Charts'
import { DataTable } from './DataTable'
import { findHistogram } from './helpers'
import { CompositionCard, EnthalpyCard, SizeCard } from './ProfileCharts'
import { ErrorNotice } from './Resource'
import { bySource, byHeavyAtoms, cumulativeAccuracy, histogram, summarize, thresholdFor, tierCounts } from './stats'
import {
  Callout, ChartCard, Explainer, Glossary, MetricCard, PageHero, ProportionBar, RankedBars, SectionHeading,
} from './StoryVisuals'
import type { Segment } from './StoryVisuals'
import { formatNumber, formatPercent, formatSigned, series, statusColor } from './theme'
import { useApi } from './useApi'

const PAGE_SIZE = 25

const FLAG_INFO: Record<string, { label: string; text: string }> = {
  validation_only: { label: 'Só validação', text: 'Valor vindo do NIST: a licença não permite usá-lo como alvo de treino de IA, apenas para conferir resultados.' },
  conflict_flag: { label: 'Conflito entre fontes', text: 'Outra fonte difere mais que max(2 × incerteza combinada, 1 kcal/mol).' },
  low_priority_only: { label: 'Só fontes de baixa prioridade', text: 'Todos os valores vêm de Stewart (usado para parametrizar o PM7) ou Yaws (mistura medidas e estimativas).' },
  stewart_reference: { label: 'Referência Stewart', text: 'Molécula usada na parametrização do PM6/PM7: o erro do PM7 nela não é um teste independente.' },
}

const TIER_COLORS = [statusColor.good, statusColor.warning, statusColor.serious, statusColor.critical]
const TIER_ICONS = ['✓', '≈', '!', '✕']

const EXP_GLOSSARY = [
  { term: 'ΔHf° (entalpia de formação)', definition: 'Energia liberada ou absorvida ao formar a molécula a partir dos elementos puros. Valores negativos indicam moléculas mais estáveis.', example: 'Metano: −17,9 kcal/mol' },
  { term: 'kcal/mol', definition: 'Unidade de energia por quantidade de moléculas. 1 kcal/mol ≈ 4,184 kJ/mol.' },
  { term: 'CBS', definition: 'Cálculo quântico de alta precisão (Complete Basis Set). Caro, mas próximo do experimento.' },
  { term: 'Precisão química', definition: 'Erro de até 1 kcal/mol — a meta clássica para previsões de energia.' },
  { term: 'Erro absoluto médio (MAE)', definition: 'Média das distâncias entre cálculo e medida, sem olhar o sinal.' },
  { term: 'Viés (erro médio)', definition: 'Média com sinal: positivo significa que o cálculo tende a dar valores acima da medida.' },
  { term: 'RMSE', definition: 'Raiz do erro quadrático médio. Pune mais os erros grandes que o MAE.' },
  { term: 'SMILES', definition: 'Forma de escrever uma molécula como texto.', example: 'CCO = etanol' },
]

function category(point: ComparisonPoint): number {
  return point.validation_only ? 1 : 0
}

function ComparisonSection({ comparison }: { comparison: ExperimentalComparison }) {
  const [query, setQuery] = useState('')
  const [source, setSource] = useState('all')
  const [validationOnly, setValidationOnly] = useState(false)
  const [conflictOnly, setConflictOnly] = useState(false)
  const [page, setPage] = useState(0)
  const sources = useMemo(() => [...new Set(comparison.points.map((point) => point.source))].sort(), [comparison.points])
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
  const stats = summarize(filtered)
  const tiers = tierCounts(filtered)
  const maxPage = Math.max(0, Math.ceil(filtered.length / PAGE_SIZE) - 1)
  const current = Math.min(page, maxPage)
  const visible = filtered.slice(current * PAGE_SIZE, (current + 1) * PAGE_SIZE)
  const resetPage = () => setPage(0)
  const p80 = thresholdFor(filtered, 80)
  const differences = histogram(filtered.map((point) => point.cbs_minus_exp), 1)
  const parity: ParityDatum[] = filtered.map((point) => ({
    x: point.H298_exp, y: point.H298_cbs, id: point.exp_id, group: category(point),
    detail: `${point.smiles} · ${point.source}`, marker: point.conflict_flag ? 'triangle' : 'circle',
  }))
  const tierSegments: Segment[] = tiers.map(({ tier, count }, index) => ({
    key: tier.key, label: `${tier.label} · ${tier.range}`, count, color: TIER_COLORS[index],
    description: `${tier.description} (kcal/mol)`, icon: TIER_ICONS[index],
  }))

  return (
    <section className="panel" aria-label="Comparação CBS × experimental">
      <SectionHeading eyebrow="Etapa 02 · comparação" title="O cálculo CBS acerta a medida?">
        Um <strong>par</strong> é uma molécula que tem valor experimental e valor CBS. A diferença entre os dois ({comparison.units.split(';')[0]}) mostra o quanto a referência calculada se afasta do laboratório.
      </SectionHeading>
      <Callout label="Como ler">
        <strong>Quanto mais perto de zero, melhor.</strong> Diferença positiva significa CBS acima do experimento. Ela não diz quem está errado: medidas também têm incerteza.
      </Callout>

      <div className="filters" role="group" aria-label="Filtros da comparação">
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
        <p className="filters__note">Indicadores, gráficos e tabela respondem aos mesmos filtros.</p>
      </div>

      <div className="metric-grid metric-grid--five" aria-label="Estatísticas dos pares filtrados">
        <MetricCard label="Pares no recorte" value={stats.n.toLocaleString('pt-BR')} note="Moléculas com CBS e experimento" accent icon="⚖" />
        <MetricCard label="Erro típico (MAE)" value={formatNumber(stats.mae)} unit="kcal/mol" note="Distância média, sem sinal" icon="↔" hint="Quanto o CBS erra, em média, para cada molécula." />
        <MetricCard label="Viés" value={formatSigned(stats.bias)} unit="kcal/mol" note="Erro médio com sinal" icon="⇡" hint={stats.bias !== null && Math.abs(stats.bias) > 1 ? `O CBS tende a ficar ${stats.bias > 0 ? 'acima' : 'abaixo'} da medida.` : 'Sem tendência forte para cima ou para baixo.'} />
        <MetricCard label="RMSE" value={formatNumber(stats.rmse)} unit="kcal/mol" note="Pesa mais os erros grandes" icon="◠" />
        <MetricCard label="Precisão química" value={formatPercent(stats.withinChemical)} note="Pares com erro ≤ 1 kcal/mol" icon="◎" tone="good" />
      </div>

      {filtered.length === 0 ? <p className="empty-state">Nenhum par corresponde aos filtros.</p> : (
        <>
          <ChartCard
            title="Classificação da concordância" wide
            subtitle="Cada molécula vai para uma faixa conforme o tamanho da diferença."
            read={<>Verde é concordância dentro da precisão química; vermelho indica casos a investigar — podem ser erros do cálculo ou da medida.</>}
          >
            <ProportionBar segments={tierSegments} label="Faixas de concordância CBS × experimental" />
          </ChartCard>

          <div className="chart-grid">
            <ChartCard
              title="Distribuição das diferenças"
              subtitle="Quantos pares caem em cada faixa de CBS − experimental (kcal/mol)."
              read={<>A linha cheia em zero é a concordância perfeita. Se o “morro” está deslocado para a direita, o CBS tende a superestimar. Aqui o viés é <strong>{formatSigned(stats.bias)} kcal/mol</strong>.</>}
              footer={differences.below + differences.above > 0 ? <p className="chart-footnote">{(differences.below + differences.above).toLocaleString('pt-BR')} pares extremos ficam fora do eixo (veja a tabela abaixo, ordenada pelo maior erro).</p> : undefined}
            >
              <HistogramChart
                histogram={{ bins: differences.bins, unit: 'kcal/mol' }}
                color={series[1]} xLabel="CBS − experimental (kcal/mol)"
                refs={[{ value: 0, label: 'zero' }, ...(stats.bias !== null ? [{ value: stats.bias, label: 'viés' }] : [])]}
              />
            </ChartCard>
            <ChartCard
              title="Acerto acumulado"
              subtitle="Que fração das moléculas fica dentro de cada margem de erro?"
              read={<>Leia como “X% das moléculas erram no máximo Y kcal/mol”. Por exemplo, 80% dos pares erram até <strong>{formatNumber(p80, 1)} kcal/mol</strong>. Quanto mais rápido a curva sobe, melhor.</>}
            >
              <CumulativeChart data={cumulativeAccuracy(filtered)} marks={[1, 5, 10]} />
            </ChartCard>
          </div>

          <div className="chart-grid">
            <ChartCard
              title="Paridade CBS × experimental"
              subtitle="Cada ponto é uma molécula: valor medido (horizontal) contra valor CBS (vertical)."
              read="Pontos sobre a linha diagonal concordam perfeitamente; as linhas tracejadas marcam ±5 kcal/mol. Pontos fora da faixa são os casos discrepantes."
              footer={<LegendRow items={[
                { label: 'Elegível para treino', color: series[0] },
                { label: 'Só validação', color: series[1] },
                { label: 'Com conflito entre fontes', marker: 'triangle' },
                { label: 'Igualdade (y = x)', marker: 'line' },
                { label: '±5 kcal/mol', marker: 'dash' },
              ]} />}
            >
              <ParityChart
                data={parity} xName="Experimental" yName="CBS"
                groups={[0, 1].flatMap((group) => (['circle', 'triangle'] as const).map((marker) => ({
                  label: `${group}-${marker}`, color: series[group], marker,
                  match: (datum: ParityDatum) => datum.group === group && datum.marker === marker,
                })))}
              />
            </ChartCard>
            <ChartCard
              title="Erro por tamanho da molécula"
              subtitle="Barras: erro absoluto médio. Linha: erro médio com sinal."
              read="Mostra se moléculas maiores são mais difíceis. Barras crescendo para a direita indicam que o erro aumenta com o tamanho."
              footer={<LegendRow items={[{ label: 'Erro absoluto médio', color: series[0] }, { label: 'Erro médio com sinal', color: series[1], marker: 'line' }]} />}
            >
              <SizeErrorChart data={byHeavyAtoms(filtered).map((item) => ({ nheavy: item.nheavy, n: item.n, mean: item.bias, mean_abs: item.mae }))} />
            </ChartCard>
          </div>

          <ChartCard
            title="Erro por fonte experimental" wide
            subtitle="Erro absoluto médio do CBS contra cada fonte (fontes com ao menos 5 pares)."
            read="Fontes com barras curtas concordam mais com o CBS. Uma barra longa sinaliza onde vale investigar a fonte ou o cálculo; a diferença sozinha não diz qual dos dois erra."
          >
            <GroupErrorChart data={bySource(filtered).map((item) => ({ key: item.key, n: item.n, mae: item.mae, bias: item.bias }))} />
          </ChartCard>

          <div className="table-heading">
            <div><span className="eyebrow eyebrow--dark">Os pares em detalhe</span><h3>Molécula por molécula</h3><p>Ordenado pela maior diferença absoluta. Valores em kcal/mol.</p></div>
          </div>
          <div className="table-wrap">
            <table className="pair-table">
              <caption className="sr-only">Pares CBS e experimental filtrados, ordenados pela diferença absoluta</caption>
              <thead><tr><th scope="col">ID exp.</th><th scope="col">ID CBS</th><th scope="col">SMILES</th><th scope="col">Fonte</th><th scope="col">Pesados</th><th scope="col">Exp.</th><th scope="col">CBS</th><th scope="col">CBS − exp</th><th scope="col">Flags</th></tr></thead>
              <tbody>{visible.map((point) => {
                const abs = Math.abs(point.cbs_minus_exp)
                const level = abs <= 1 ? 'good' : abs <= 5 ? 'warning' : abs <= 10 ? 'serious' : 'critical'
                return (
                  <tr key={point.exp_id}>
                    <td className="id-cell">{point.exp_id}</td><td className="id-cell">{point.mol_id}</td><td className="smiles">{point.smiles}</td><td>{point.source}</td>
                    <td className="num-cell">{point.nheavy}</td><td className="num-cell">{formatNumber(point.H298_exp)}</td><td className="num-cell">{formatNumber(point.H298_cbs)}</td>
                    <td className="num-cell"><span className={`difference difference--${level}`}>{formatSigned(point.cbs_minus_exp)}</span></td>
                    <td className="flags-cell">{[point.validation_only && 'só validação', point.conflict_flag && 'conflito', point.low_priority_only && 'baixa prioridade'].filter(Boolean).join(', ') || '—'}</td>
                  </tr>
                )
              })}</tbody>
            </table>
          </div>
          <div className="pagination">
            <span>Mostrando {current * PAGE_SIZE + 1}–{Math.min((current + 1) * PAGE_SIZE, filtered.length)} de {filtered.length.toLocaleString('pt-BR')}</span>
            <button type="button" disabled={current === 0} onClick={() => setPage(current - 1)}>Anterior</button>
            <button type="button" disabled={current >= maxPage} onClick={() => setPage(current + 1)}>Próxima</button>
          </div>
        </>
      )}
    </section>
  )
}

function Distributions({ summary, profile }: { summary: ExperimentalSummary; profile: DatasetProfile | null }) {
  const familySegments: Segment[] = Object.entries(summary.by_source_family)
    .sort((a, b) => b[1] - a[1])
    .map(([name, count], index) => ({ key: name, label: name, count, color: series[index % series.length] }))
  const uncertainty = profile ? findHistogram(profile, 'uncertainty') : undefined
  const reported = summary.uncertainty_kind.reported ?? 0
  return (
    <>
      <SectionHeading eyebrow="Etapa 01 · o que há na tabela" title="Como os dados experimentais se distribuem">
        Cada gráfico responde a uma pergunta simples sobre as {summary.rows.toLocaleString('pt-BR')} moléculas medidas.
      </SectionHeading>
      {profile && (
        <div className="chart-grid">
          <EnthalpyCard
            histogram={findHistogram(profile, 'H298_exp')} title="Quão estáveis são as moléculas?"
            subtitle="Distribuição da entalpia de formação medida (kcal/mol)."
          />
          <SizeCard profile={profile} />
        </div>
      )}
      <div className="chart-grid">
        <ChartCard
          title="De onde vêm as medidas?"
          subtitle="Famílias de bancos de dados; cada molécula tem uma fonte selecionada."
          read="ATcT (Active Thermochemical Tables) é a referência mais rigorosa; Bains reúne compilações clássicas (Pedley, Yaws…); NIST é o WebBook, usado só para validação; RMG traz o conjunto de referência do projeto RMG."
        >
          <DonutChart segments={familySegments} centerValue={summary.rows.toLocaleString('pt-BR')} centerCaption="moléculas" label="Moléculas por família de fonte" />
        </ChartCard>
        <ChartCard
          title="Fontes individuais"
          subtitle="As seis fontes que mais contribuem; o restante é agrupado."
          read="Ajuda a ver se o conjunto depende demais de um único banco."
        >
          <RankedBars counts={summary.by_source} total={summary.rows} label="Moléculas por fonte" />
        </ChartCard>
      </div>
      <div className="chart-grid">
        {profile && <CompositionCard profile={profile} />}
        <ChartCard
          title="Quão confiáveis são as medidas?"
          subtitle={`${reported.toLocaleString('pt-BR')} moléculas (${formatPercent(summary.rows ? reported / summary.rows : null)}) têm incerteza informada pela fonte.`}
          read="A incerteza é o “±” da medida. Valores pequenos indicam medidas mais precisas; muitas fontes clássicas não informam incerteza."
        >
          <ProportionBar
            label="Incerteza informada"
            segments={[
              { key: 'reported', label: 'Com incerteza informada', count: reported, color: series[0] },
              { key: 'none', label: 'Sem incerteza informada', count: summary.uncertainty_kind.none ?? 0, color: 'var(--status-neutral)' },
            ]}
          />
          {uncertainty && uncertainty.n > 0 && (
            <HistogramChart histogram={uncertainty} color={series[0]} xLabel="Incerteza (± kcal/mol)" height={190}
              refs={uncertainty.median !== null ? [{ value: uncertainty.median, label: 'mediana' }] : []} />
          )}
        </ChartCard>
      </div>
      <ChartCard
        title="Observações sobre os dados" wide
        subtitle="Marcadores de qualidade. Uma molécula pode ter mais de um, então não somam 100%."
      >
        <ul className="flag-list">
          {Object.entries(summary.flags).map(([flag, count]) => (
            <li key={flag}>
              <div className="flag-list__head"><strong>{FLAG_INFO[flag]?.label ?? flag}</strong><span>{count.toLocaleString('pt-BR')} · {formatPercent(summary.rows ? count / summary.rows : null)}</span></div>
              <span className="flag-list__track"><span style={{ width: `${summary.rows ? (count / summary.rows) * 100 : 0}%` }} /></span>
              {FLAG_INFO[flag] && <p>{FLAG_INFO[flag].text}</p>}
            </li>
          ))}
        </ul>
      </ChartCard>
    </>
  )
}

export function ExperimentalView() {
  const summary = useApi<ExperimentalSummary>('/api/experimental/summary')
  const comparison = useApi<ExperimentalComparison>('/api/experimental/comparison')
  const profile = useApi<DatasetProfile>('/api/tables/exp/profile')
  const sum = summary.state.status === 'ready' ? summary.state.data : null
  const prof = profile.state.status === 'ready' ? profile.state.data : null
  const pairs = comparison.state.status === 'ready' ? comparison.state.data.points.length : null
  return (
    <div className="view">
      <PageHero
        variant="exp" eyebrow="02 · Tabela experimental" title="O que o laboratório" accent="mediu de fato."
        lead="Valores medidos de entalpia de formação: o gabarito do mundo real contra o qual comparamos os cálculos."
        onRefresh={() => { summary.reload(); comparison.reload(); profile.reload() }}
      />
      <Explainer
        title="O que é a tabela experimental?"
        intro="Uma coleção de entalpias de formação medidas em laboratório, em fase gasosa a 25 °C, para moléculas feitas apenas de C, H, O e N. Ela é montada a partir de vários bancos públicos e cada molécula recebe um único valor selecionado."
        analogy="Se o CBS é uma calculadora muito boa, a tabela experimental é a balança calibrada que diz se a calculadora está certa."
        facts={[
          { icon: '🧪', title: 'Origem', text: 'ATcT, NIST WebBook, compilações de Bains et al. (2022) e o conjunto de referência do RMG.' },
          { icon: '🎯', title: 'Para que serve', text: 'Validar os cálculos e, quando elegível, servir de alvo de treino.' },
          { icon: '⚠️', title: 'Cuidados', text: 'Valores do NIST são “só validação” (licença) e nunca entram no treino; conflitos entre fontes ficam sinalizados.' },
        ]}
      />
      {summary.state.status === 'loading' && <p role="status" className="loading-state">Carregando fontes…</p>}
      {summary.state.status === 'error' && <ErrorNotice message={summary.state.message} retry={summary.reload} />}
      {profile.state.status === 'error' && <ErrorNotice message={profile.state.message} retry={profile.reload} />}
      {sum && (sum.rows === 0 ? <p className="empty-state">Nenhum valor experimental disponível.</p> : (
        <>
          <div className="metric-grid metric-grid--four">
            <div className="metric-card metric-card--accent">
              <p className="metric-card__label">Moléculas medidas</p>
              <p className="metric metric-card__value">{sum.rows.toLocaleString('pt-BR')}</p>
              <p className="metric-card__note">ΔHf° gasoso a 298,15 K</p>
            </div>
            <MetricCard label="Fontes distintas" value={Object.keys(sum.by_source).length.toLocaleString('pt-BR')} note={`${Object.keys(sum.by_source_family).length} famílias de bancos`} icon="📚" />
            <MetricCard label="Com incerteza" value={formatPercent(sum.rows ? (sum.uncertainty_kind.reported ?? 0) / sum.rows : null)} note="Fonte informa o ± da medida" icon="±" />
            <MetricCard label="Pares com CBS" value={pairs === null ? '—' : pairs.toLocaleString('pt-BR')} note="Moléculas comparáveis ao cálculo" icon="⚖" />
          </div>
          <Distributions summary={sum} profile={prof} />
        </>
      ))}
      {comparison.state.status === 'loading' && <p role="status" className="loading-state">Carregando comparação CBS…</p>}
      {comparison.state.status === 'error' && <section className="panel" aria-label="Comparação CBS"><ErrorNotice message={comparison.state.message} retry={comparison.reload} /></section>}
      {comparison.state.status === 'ready' && (comparison.state.data.points.length === 0
        ? <section className="panel" aria-label="Comparação CBS"><p className="empty-state">Nenhum par CBS × experimental disponível.</p></section>
        : <ComparisonSection comparison={comparison.state.data} />)}
      {sum && sum.rows > 0 && (
        <section className="panel" aria-label="Tabela experimental completa">
          <SectionHeading eyebrow="Etapa 03 · consulta" title="A tabela completa">Pesquise por nome, fórmula, SMILES ou fonte.</SectionHeading>
          <DataTable table="exp" title="Valores experimentais" searchHint="Ex.: ethanol, C2H6O, atct" defaultSort={{ column: 'exp_id', desc: false }} />
        </section>
      )}
      <Glossary terms={EXP_GLOSSARY} />
    </div>
  )
}

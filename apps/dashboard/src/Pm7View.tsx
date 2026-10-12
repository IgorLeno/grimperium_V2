import type { DatasetProfile } from './api'
import { MissingTable } from './CbsView'
import { LegendRow, ParityChart, SizeErrorChart } from './Charts'
import { DataTable } from './DataTable'
import { findHistogram } from './helpers'
import { CompositionCard, EnthalpyCard, SizeCard } from './ProfileCharts'
import {
  Callout, ChartCard, DeltaExplainer, Explainer, Glossary, MetricCard, PageHero, ProportionBar, RingGauge, SectionHeading,
} from './StoryVisuals'
import { statusSegments } from './helpers'
import type { Segment } from './StoryVisuals'
import { formatNumber, formatSigned, series, statusColor } from './theme'
import { useApi } from './useApi'

function outcomeSegments(counts: Record<string, number> | undefined): Segment[] {
  const palette: Record<string, string> = {
    SUCCESS: statusColor.good, TIMEOUT: statusColor.warning, ERROR: statusColor.critical, NOT_ATTEMPTED: statusColor.neutral,
  }
  return Object.entries(counts ?? {}).sort((a, b) => b[1] - a[1]).map(([key, count]) => ({
    key, label: key === 'unset' ? 'sem registro' : key, count, color: palette[key] ?? 'var(--status-neutral)',
  }))
}

const PM7_GLOSSARY = [
  { term: 'PM7', definition: 'Método semiempírico: usa parâmetros ajustados para calcular energias em segundos. Rápido, mas com erros de vários kcal/mol.' },
  { term: 'CREST', definition: 'Programa que procura as conformações de menor energia antes do cálculo final.' },
  { term: 'MOPAC', definition: 'Programa que executa o cálculo PM7.' },
  { term: 'Δ (delta)', definition: 'Diferença CBS − PM7. É o que o modelo de IA aprende a prever.' },
  { term: 'Flag de qualidade CBS', definition: 'Marca referências CBS suspeitas (SUSPECT) para que não contaminem o treino.' },
]

export function Pm7View() {
  const { state, reload } = useApi<DatasetProfile>('/api/tables/pm7/profile')
  const profile = state.status === 'ready' ? state.data : null
  const status = profile?.categories.status ?? {}
  const done = status.OK ?? 0
  const delta = profile ? findHistogram(profile, 'delta_cbs_pm7') : undefined
  const crestTime = profile ? findHistogram(profile, 'crest_time_s') : undefined
  const mopacTime = profile ? findHistogram(profile, 'mopac_time_s') : undefined
  const quality = profile?.categories.cbs_quality_flag ?? {}
  return (
    <div className="view">
      <PageHero
        variant="pm7" eyebrow="04 · Tabela PM7" title="O cálculo rápido" accent="e o erro que ele comete."
        lead="Para cada molécula da referência CBS, o pipeline CREST + MOPAC calcula a entalpia com PM7. A diferença para o CBS é o alvo do modelo."
        onRefresh={reload}
      />
      <Explainer
        title="O que é a tabela PM7?"
        intro="É a tabela construída pelo pipeline legado crest_pm7: uma linha por molécula, com estado do processamento, resultado do CREST e do MOPAC, entalpia PM7, descritores e a diferença para o CBS. O protocolo canônico agora é o Semi-Imperium; esta tabela é mostrada como legado."
        analogy="É o boletim do aluno rápido: para cada questão (molécula), a resposta dele (PM7) e a distância até o gabarito (CBS)."
        facts={[
          { icon: '⚙️', title: 'Como é gerada', text: 'CREST busca conformações → MOPAC calcula PM7 → o resultado é gravado com o status da molécula.' },
          { icon: '🎯', title: 'Para que serve', text: 'Fornecer Δ = CBS − PM7, o alvo do treino, e descritores para o modelo.' },
          { icon: '🧹', title: 'Filtro de treino', text: 'Só entram linhas com status OK e qualidade CBS OK.' },
        ]}
      />
      <section className="panel" aria-label="Como a correção funciona">
        <SectionHeading eyebrow="Ideia central" title="Por que calcular PM7 se temos CBS?">
          Porque o CBS não escala para moléculas novas. O PM7 sim — e o modelo corrige o erro dele.
        </SectionHeading>
        <DeltaExplainer />
      </section>
      {state.status === 'loading' && <p role="status" className="loading-state">Carregando tabela PM7…</p>}
      {state.status === 'error' && <MissingTable message={state.message} retry={reload} file="data/thermo_pm7.csv" />}
      {profile && (
        <>
          <section className="panel" aria-label="Andamento do pipeline">
            <SectionHeading eyebrow="Andamento" title="Quanto já foi calculado?" badge="Protocolo legado: crest_pm7">
              Cada molécula passa por estados até ficar OK. A barra mostra a fila inteira.
            </SectionHeading>
            <div className="progress-layout">
              <RingGauge share={profile.rows ? done / profile.rows : null} label="Moléculas concluídas" caption={`${done.toLocaleString('pt-BR')} de ${profile.rows.toLocaleString('pt-BR')}`} color={statusColor.good} />
              <ProportionBar segments={statusSegments(status)} total={profile.rows} label="Estados das moléculas" />
            </div>
          </section>

          <div className="metric-grid metric-grid--four">
            <MetricCard label="Moléculas com Δ" value={(delta?.n ?? 0).toLocaleString('pt-BR')} note="Têm CBS e PM7 calculados" accent icon="Δ" />
            <MetricCard label="Erro típico do PM7" value={formatNumber(delta?.median != null ? Math.abs(delta.median) : null, 1)} unit="kcal/mol" note="|mediana| de CBS − PM7" icon="↔" hint="Quanto o modelo precisa corrigir numa molécula típica." />
            <MetricCard label="Qualidade CBS OK" value={(quality.OK ?? 0).toLocaleString('pt-BR')} note={`${(quality.SUSPECT ?? 0).toLocaleString('pt-BR')} marcadas como suspeitas`} icon="🛡" />
            <MetricCard label="Tempo típico" value={formatNumber(crestTime?.median ?? null, 0)} unit="s" note={`CREST · MOPAC ${formatNumber(mopacTime?.median ?? null, 0)} s (medianas)`} icon="⏱" />
          </div>

          <SectionHeading eyebrow="O erro do PM7" title="Onde o cálculo rápido erra">
            O erro não é aleatório: tem tamanho e direção típicos — por isso dá para aprendê-lo.
          </SectionHeading>
          <div className="chart-grid">
            <EnthalpyCard
              histogram={delta} title="Distribuição do erro (Δ = CBS − PM7)" color={series[1]}
              subtitle="Quantas moléculas têm cada tamanho de correção."
              refs={delta?.mean != null ? [{ value: 0, label: 'zero' }, { value: delta.mean, label: 'média' }] : []}
              read={<>Se o PM7 fosse perfeito, todas as barras estariam em zero. O deslocamento médio é <strong>{formatSigned(delta?.mean ?? null)} kcal/mol</strong>: uma parte do erro é sistemática, e o modelo aprende o resto a partir da estrutura.</>}
            />
            <ChartCard
              title="Paridade CBS × PM7"
              subtitle={`Cada ponto é uma molécula${profile.pairs_total > profile.pairs.length ? ` (amostra de ${profile.pairs.length.toLocaleString('pt-BR')} de ${profile.pairs_total.toLocaleString('pt-BR')})` : ''}.`}
              read="Pontos na diagonal: PM7 igual ao CBS. A distância vertical até a diagonal é o Δ que o modelo vai prever."
              footer={<LegendRow items={[{ label: 'Molécula', color: series[0] }, { label: 'Igualdade (y = x)', marker: 'line' }, { label: '±5 kcal/mol', marker: 'dash' }]} />}
            >
              {profile.pairs.length ? (
                <ParityChart
                  data={profile.pairs.map((pair) => ({ x: pair.x, y: pair.y, id: pair.id, group: 0, detail: pair.nheavy !== null ? `${pair.nheavy} átomos pesados` : undefined }))}
                  groups={[{ label: 'Moléculas', color: series[0], match: () => true }]} xName="CBS" yName="PM7"
                />
              ) : <p className="empty-state">Ainda não há moléculas com CBS e PM7 juntos.</p>}
            </ChartCard>
          </div>
          <div className="chart-grid">
            <ChartCard
              title="O erro cresce com o tamanho?"
              subtitle="Barras: erro absoluto médio por número de átomos pesados. Linha: média com sinal."
              read="Se as barras crescem para a direita, moléculas maiores acumulam mais erro do PM7 — um motivo para o tamanho ser uma das variáveis do modelo."
              footer={<LegendRow items={[{ label: 'Erro absoluto médio', color: series[0] }, { label: 'Erro médio com sinal', color: series[1], marker: 'line' }]} />}
            >
              {profile.by_nheavy.length ? <SizeErrorChart data={profile.by_nheavy} /> : <p className="empty-state">Sem moléculas com Δ e tamanho.</p>}
            </ChartCard>
            <ChartCard title="Resultado de cada programa" subtitle="Como terminaram as execuções do CREST e do MOPAC.">
              <h4 className="mini-heading">CREST (busca conformacional)</h4>
              <ProportionBar segments={outcomeSegments(profile.categories.crest_status)} label="Resultado do CREST" />
              <h4 className="mini-heading">MOPAC (cálculo PM7)</h4>
              <ProportionBar segments={outcomeSegments(profile.categories.mopac_status)} label="Resultado do MOPAC" />
              <h4 className="mini-heading">Qualidade da referência CBS</h4>
              <ProportionBar segments={Object.entries(quality).map(([key, count]) => ({
                key, label: key === 'unset' ? 'sem flag' : key, count,
                color: key === 'OK' ? statusColor.good : key === 'SUSPECT' ? statusColor.warning : statusColor.neutral,
              }))} label="Flag de qualidade CBS" />
            </ChartCard>
          </div>
          <div className="chart-grid">
            <SizeCard profile={profile} color={series[2]} />
            <CompositionCard profile={profile} />
          </div>
          {(crestTime?.n ?? 0) + (mopacTime?.n ?? 0) > 0 && (
            <div className="chart-grid">
              <EnthalpyCard histogram={crestTime} title="Quanto tempo leva o CREST?" subtitle="Tempo de busca conformacional por molécula (s)." color={series[2]} />
              <EnthalpyCard histogram={mopacTime} title="Quanto tempo leva o MOPAC?" subtitle="Tempo do cálculo PM7 por molécula (s)." color={series[3]} />
            </div>
          )}
          {profile.missing_columns.length > 0 && (
            <Callout label="Colunas ausentes" tone="warning">Este CSV não tem: {profile.missing_columns.join(', ')}. Os gráficos correspondentes ficam vazios.</Callout>
          )}
          <section className="panel" aria-label="Tabela PM7 completa">
            <SectionHeading eyebrow="Consulta" title="A tabela completa">Pesquise por ID, SMILES ou estado; ordene por Δ para ver os maiores erros.</SectionHeading>
            <DataTable table="pm7" title="Resultados CREST + PM7" searchHint="Ex.: OK, Pending, CCO" defaultSort={{ column: 'mol_id', desc: false }} />
          </section>
        </>
      )}
      <Glossary terms={PM7_GLOSSARY} />
    </div>
  )
}

import type { ReactNode } from 'react'
import type { DatasetProfile, DatasetProgress, ExperimentalComparison } from './api'
import { statusSegments } from './helpers'
import { summarize } from './stats'
import {
  DeltaExplainer, Glossary, MetricCard, PageHero, PipelineDiagram, ProportionBar, RingGauge, SectionHeading,
} from './StoryVisuals'
import { formatNumber, formatPercent, formatSigned, series, statusColor } from './theme'
import { useApi } from './useApi'

export type ViewKey = 'summary' | 'experimental' | 'cbs' | 'pm7' | 'training' | 'learning' | 'calc'

const CORE_GLOSSARY = [
  { term: 'Entalpia de formação (ΔHf°)', definition: 'A energia de “montar” a molécula a partir dos elementos. É a propriedade que o Grimperium prevê.' },
  { term: 'PM7', definition: 'Cálculo rápido e aproximado (segundos por molécula).' },
  { term: 'CBS', definition: 'Cálculo lento e muito preciso; serve de gabarito.' },
  { term: 'Experimental', definition: 'Valor medido em laboratório; confere se os cálculos batem com a realidade.' },
  { term: 'Δ-learning', definition: 'Estratégia em que a IA aprende só o erro do PM7 (Δ = CBS − PM7), e não a energia inteira. Precisa de menos dados.' },
  { term: 'Semi-Imperium', definition: 'O novo pipeline canônico de cálculo (RDKit → CREST → seleção → MOPAC), com registro completo de proveniência.' },
]

function n(value: number | null): string {
  return value === null ? '—' : value.toLocaleString('pt-BR')
}

const STORE_STATES: Record<string, { label: string; color: string; icon: string }> = {
  verified: { label: 'Verificado', color: statusColor.good, icon: '✓' },
  unverified: { label: 'Sem verificação', color: series[0], icon: '○' },
  running: { label: 'Rodando', color: statusColor.info, icon: '▶' },
  pending: { label: 'Pendente', color: statusColor.neutral, icon: '…' },
  saddle: { label: 'Ponto de sela', color: statusColor.warning, icon: '!' },
  failed: { label: 'Falhou', color: statusColor.critical, icon: '✕' },
}

function storeSegments(states: Record<string, number>) {
  return Object.entries(states).filter(([, value]) => value > 0).map(([key, value]) => ({
    key, label: STORE_STATES[key]?.label ?? key, icon: STORE_STATES[key]?.icon, count: value, color: STORE_STATES[key]?.color ?? statusColor.neutral,
  }))
}

function TableNode({ title, value, caption, tone, onOpen, unavailable }: {
  title: string
  value: string
  caption: string
  tone: 'exp' | 'cbs' | 'pm7'
  onOpen: () => void
  unavailable: boolean
}) {
  return (
    <button type="button" className={`data-node data-node--${tone}`} onClick={onOpen}>
      <span className="data-node__title">{title}</span>
      <strong>{value}</strong>
      <span className="data-node__caption">{unavailable ? 'arquivo não encontrado neste computador' : caption}</span>
      <span className="data-node__open">Abrir tabela →</span>
    </button>
  )
}

function Link({ children }: { children: ReactNode }) {
  return <span className="data-link"><span aria-hidden="true" className="data-link__line" /><span className="data-link__label">{children}</span></span>
}

export function SummaryView({ onNavigate }: { onNavigate: (view: ViewKey) => void }) {
  const progress = useApi<DatasetProgress>('/api/dataset/progress')
  const exp = useApi<DatasetProfile>('/api/tables/exp/profile')
  const cbs = useApi<DatasetProfile>('/api/tables/cbs/profile')
  const comparison = useApi<ExperimentalComparison>('/api/experimental/comparison')
  const pm7 = useApi<DatasetProfile>('/api/tables/pm7/profile')

  const expRows = exp.state.status === 'ready' ? exp.state.data.rows : null
  const cbsRows = cbs.state.status === 'ready' ? cbs.state.data.rows : null
  const legacy = progress.state.status === 'ready' ? progress.state.data.legacy : null
  const store = progress.state.status === 'ready' ? progress.state.data.semi_imperium : null
  const pm7Profile = pm7.state.status === 'ready' ? pm7.state.data : null
  const pm7Total = legacy?.total ?? pm7Profile?.rows ?? null
  const pm7Done = legacy?.status.OK ?? pm7Profile?.categories.status?.OK ?? null
  const stats = comparison.state.status === 'ready' ? summarize(comparison.state.data.points) : null
  const delta = pm7Profile?.histograms.find((item) => item.column === 'delta_cbs_pm7')

  const reloadAll = () => { progress.reload(); exp.reload(); cbs.reload(); comparison.reload(); pm7.reload() }

  return (
    <div className="view">
      <PageHero
        eyebrow="Resumo do módulo Delta Entalpia" title="De moléculas a" accent="conhecimento confiável."
        lead="O Grimperium prevê a energia de formação de moléculas combinando um cálculo rápido (PM7) com uma correção aprendida por IA a partir de cálculos de alta precisão (CBS)."
        onRefresh={reloadAll}
      />

      <section className="panel summary-text" aria-label="Resumo em texto">
        <SectionHeading eyebrow="Em poucas palavras" title="Onde o projeto está hoje" />
        <div className="summary-text__body">
          <p>
            A <strong>referência CBS</strong> {cbsRows === null ? 'não está disponível neste computador' : <>reúne <strong>{n(cbsRows)}</strong> geometrias de moléculas CHON calculadas com alta precisão</>}.
            {' '}A <strong>tabela experimental</strong> {expRows === null ? 'não pôde ser lida' : <>tem <strong>{n(expRows)}</strong> moléculas medidas em laboratório</>}
            {stats && <>, das quais <strong>{n(stats.n)}</strong> também têm valor CBS</>}.
          </p>
          {stats && (
            <p>
              Nesses pares, o CBS erra em média <strong>{formatNumber(stats.mae, 1)} kcal/mol</strong>, com tendência de ficar
              {' '}{stats.bias !== null && stats.bias >= 0 ? 'acima' : 'abaixo'} da medida ({formatSigned(stats.bias, 1)} kcal/mol);
              {' '}<strong>{formatPercent(stats.withinChemical)}</strong> dos pares atingem a precisão química (≤ 1 kcal/mol).
            </p>
          )}
          <p>
            O pipeline legado <strong>CREST + PM7</strong> {pm7Total === null ? 'não tem CSV disponível neste computador' : <>concluiu <strong>{n(pm7Done)}</strong> de <strong>{n(pm7Total)}</strong> moléculas ({formatPercent(pm7Total ? (pm7Done ?? 0) / pm7Total : null)})</>}
            {delta && delta.n > 0 && <>; nelas, a correção típica que o modelo precisa aprender é de {formatNumber(delta.median !== null ? Math.abs(delta.median) : null, 1)} kcal/mol</>}.
            {' '}O novo pipeline <strong>Semi-Imperium</strong> {store === null ? 'não pôde ser consultado' : store.exists ? <>registra <strong>{n(store.molecules)}</strong> moléculas e <strong>{n(store.calculations)}</strong> cálculos</> : 'ainda não tem execuções registradas'}.
          </p>
        </div>
      </section>

      <section className="panel" aria-label="Como funciona">
        <SectionHeading eyebrow="Como funciona" title="O caminho de uma molécula">Cinco passos, do texto da estrutura até a previsão corrigida.</SectionHeading>
        <PipelineDiagram />
      </section>

      <section className="panel" aria-label="Mapa dos dados">
        <SectionHeading eyebrow="Mapa dos dados" title="Três tabelas, uma pergunta">Cada tabela tem um papel. Clique para explorar.</SectionHeading>
        <div className="data-map">
          <TableNode title="Experimental" value={n(expRows)} caption="moléculas medidas" tone="exp" onOpen={() => onNavigate('experimental')} unavailable={exp.state.status === 'error'} />
          <Link>{stats ? `${n(stats.n)} pares comparados` : 'comparação'}<br /><small>valida o CBS</small></Link>
          <TableNode title="Referência CBS" value={n(cbsRows)} caption="geometrias de alta precisão" tone="cbs" onOpen={() => onNavigate('cbs')} unavailable={cbs.state.status === 'error'} />
          <Link>Δ = CBS − PM7<br /><small>alvo do modelo</small></Link>
          <TableNode title="PM7 (legado)" value={n(pm7Total)} caption="moléculas na fila do pipeline" tone="pm7" onOpen={() => onNavigate('pm7')} unavailable={pm7Total === null && pm7.state.status === 'error'} />
        </div>
      </section>

      <div className="metric-grid metric-grid--four">
        <MetricCard label="Medidas experimentais" value={n(expRows)} note="Gabarito do mundo real" icon="🧪" accent />
        <MetricCard label="Referência CBS" value={n(cbsRows)} note="Gabarito calculado" icon="🏛" />
        <MetricCard label="Erro CBS × experimento" value={formatNumber(stats?.mae ?? null, 1)} unit="kcal/mol" note="Erro absoluto médio" icon="⚖" />
        <MetricCard label="PM7 concluído" value={formatPercent(pm7Total ? (pm7Done ?? 0) / pm7Total : null)} note={`${n(pm7Done)} de ${n(pm7Total)} moléculas`} icon="⚡" />
      </div>

      <div className="chart-grid">
        <section className="panel" aria-label="Dataset PM7">
          <SectionHeading eyebrow="Pipeline legado" title="Fila do CREST + PM7" badge={legacy ? `Protocolo legado: ${legacy.protocol}` : undefined} />
          {progress.state.status === 'loading' && <p role="status" className="loading-state">Carregando progresso…</p>}
          {progress.state.status === 'error' && <p className="empty-state" role="alert">{progress.state.message}</p>}
          {legacy && (legacy.total === 0 ? <p className="empty-state">Nenhuma linha no dataset.</p> : (
            <div className="progress-layout progress-layout--stack">
              <RingGauge share={legacy.total ? (legacy.status.OK ?? 0) / legacy.total : null} label="Concluídas" caption={`${(legacy.status.OK ?? 0).toLocaleString('pt-BR')} de ${legacy.total.toLocaleString('pt-BR')} linhas`} color={statusColor.good} />
              <ProportionBar segments={statusSegments(legacy.status)} total={legacy.total} label="Estados do dataset PM7" />
              <h4 className="mini-heading">Qualidade da referência CBS</h4>
              <p className="subtle">Flag de qualidade, independente do estado do cálculo</p>
              <ProportionBar label="Flag de qualidade CBS" segments={Object.entries(legacy.cbs_quality_flag).map(([key, value]) => ({
                key, label: key === 'unset' ? 'sem flag' : key, count: value,
                color: key === 'OK' ? statusColor.good : key === 'SUSPECT' ? statusColor.warning : statusColor.neutral,
              }))} />
            </div>
          ))}
        </section>
        <section className="panel" aria-label="Semi-Imperium">
          <SectionHeading eyebrow="Pipeline canônico" title="Store Semi-Imperium">Contagens separadas do CSV legado.</SectionHeading>
          {store && (!store.exists ? <p className="empty-state">Store ainda não existe. A primeira rotina executada cria o registro.</p> : (
            <>
              <div className="mini-metrics">
                <div><strong>{n(store.molecules)}</strong><span>moléculas</span></div>
                <div><strong>{n(store.runs)}</strong><span>execuções</span></div>
                <div><strong>{n(store.calculations)}</strong><span>cálculos</span></div>
              </div>
              {store.calculations === 0 ? <p className="empty-state">Nenhum cálculo registrado.</p> : Object.values(store.states).every((value) => value === 0) ? (
                <p className="empty-state">Há {n(store.calculations)} cálculo(s) gravado(s), mas nenhum ligado a uma execução ({n(store.runs)} execuções) — por isso não há “estado mais recente” para mostrar.</p>
              ) : (
                <ProportionBar label="Estado mais recente dos cálculos" segments={storeSegments(store.states)} />
              )}
            </>
          ))}
          <button type="button" className="ghost-button" onClick={() => onNavigate('training')}>Acompanhar o treinamento →</button>
        </section>
      </div>

      <section className="panel" aria-label="A ideia do Δ-learning">
        <SectionHeading eyebrow="A ideia central" title="Rápido + correção = preciso" />
        <DeltaExplainer />
      </section>

      <section className="panel roadmap" aria-label="Próximos passos">
        <SectionHeading eyebrow="Em construção" title="O que vem a seguir" />
        <ol>
          <li><strong>Treinamento em lote.</strong> A tela Treinamento já aciona e acompanha a rotina; o lote completo das moléculas pendentes ainda depende da validação do protocolo 2.</li>
          <li><strong>Planilha online.</strong> Espelho somente leitura das tabelas no Google Sheets, atualizado a partir dos CSVs locais.</li>
          <li><strong>Treino e avaliação.</strong> Painel do modelo Δ-learning com o quality gate (MAE ≤ 3,5, R² ≥ 0,97, RMSE ≤ 5,0 kcal/mol).</li>
        </ol>
      </section>

      <Glossary terms={CORE_GLOSSARY} title="Os seis termos essenciais" />
    </div>
  )
}

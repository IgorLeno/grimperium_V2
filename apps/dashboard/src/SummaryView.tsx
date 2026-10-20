import type { ReactNode } from 'react'
import type { DatasetProfile, DatasetProgress, ExperimentalComparison } from './api'
import { queueSegments } from './helpers'
import { summarize } from './stats'
import {
  DeltaExplainer, Glossary, MetricCard, PageHero, PipelineDiagram, ProportionBar, RingGauge, SectionHeading,
} from './StoryVisuals'
import { formatNumber, formatPercent, formatSigned, series, statusColor } from './theme'
import { useApi } from './useApi'

export type ViewKey = 'summary' | 'experimental' | 'cbs' | 'training' | 'learning' | 'calc'

const CORE_GLOSSARY = [
  { term: 'Entalpia de formação (ΔHf°)', definition: 'A energia de “montar” a molécula a partir dos elementos. É a propriedade que o Grimperium prevê.' },
  { term: 'PM7', definition: 'Cálculo rápido e aproximado (segundos por molécula).' },
  { term: 'CBS', definition: 'Cálculo lento e muito preciso; serve de gabarito.' },
  { term: 'Experimental', definition: 'Valor medido em laboratório; confere se os cálculos batem com a realidade.' },
  { term: 'Δ-learning', definition: 'Estratégia em que a IA aprende só o erro do PM7 (Δ = CBS − PM7), e não a energia inteira. Precisa de menos dados.' },
  { term: 'Semi-Imperium', definition: 'O pipeline canônico de cálculo (RDKit → CREST → seleção → MOPAC), com registro completo de proveniência.' },
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
  tone: 'exp' | 'cbs'
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

  const expRows = exp.state.status === 'ready' ? exp.state.data.rows : null
  const cbsRows = cbs.state.status === 'ready' ? cbs.state.data.rows : null
  const queue = progress.state.status === 'ready' ? progress.state.data.queue : null
  const store = progress.state.status === 'ready' ? progress.state.data.semi_imperium : null
  const stats = comparison.state.status === 'ready' ? summarize(comparison.state.data.points) : null

  const reloadAll = () => { progress.reload(); exp.reload(); cbs.reload(); comparison.reload() }

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
            O pipeline <strong>Semi-Imperium</strong> (protocolo 2: RDKit → CREST → seleção → MOPAC PM7) {queue === null ? 'não pôde ser consultado' : <>concluiu <strong>{n(queue.ok)}</strong> de <strong>{n(queue.total)}</strong> moléculas da fila ({formatPercent(queue.total ? queue.ok / queue.total : null)})</>}
            {store === null ? '' : store.exists ? <>; o store registra <strong>{n(store.molecules)}</strong> moléculas e <strong>{n(store.calculations)}</strong> cálculos</> : '; ainda não há execuções registradas no store'}.
          </p>
        </div>
      </section>

      <section className="panel" aria-label="Como funciona">
        <SectionHeading eyebrow="Como funciona" title="O caminho de uma molécula">Cinco passos, do texto da estrutura até a previsão corrigida.</SectionHeading>
        <PipelineDiagram />
      </section>

      <section className="panel" aria-label="Mapa dos dados">
        <SectionHeading eyebrow="Mapa dos dados" title="Duas tabelas, uma pergunta">Cada tabela tem um papel. Clique para explorar.</SectionHeading>
        <div className="data-map">
          <TableNode title="Experimental" value={n(expRows)} caption="moléculas medidas" tone="exp" onOpen={() => onNavigate('experimental')} unavailable={exp.state.status === 'error'} />
          <Link>{stats ? `${n(stats.n)} pares comparados` : 'comparação'}<br /><small>valida o CBS</small></Link>
          <TableNode title="Referência CBS" value={n(cbsRows)} caption="geometrias de alta precisão" tone="cbs" onOpen={() => onNavigate('cbs')} unavailable={cbs.state.status === 'error'} />
        </div>
      </section>

      <div className="metric-grid">
        <MetricCard label="Medidas experimentais" value={n(expRows)} note="Gabarito do mundo real" icon="🧪" accent />
        <MetricCard label="Referência CBS" value={n(cbsRows)} note="Gabarito calculado" icon="🏛" />
        <MetricCard label="Erro CBS × experimento" value={formatNumber(stats?.mae ?? null, 1)} unit="kcal/mol" note="Erro absoluto médio" icon="⚖" />
      </div>

      <div className="chart-grid">
        <section className="panel" aria-label="Fila do treinamento">
          <SectionHeading eyebrow="Pipeline canônico" title="Fila do treinamento">Uma molécula por SMILES, tirada da referência CBS.</SectionHeading>
          {progress.state.status === 'loading' && <p role="status" className="loading-state">Carregando progresso…</p>}
          {progress.state.status === 'error' && <p className="empty-state" role="alert">{progress.state.message}</p>}
          {queue && (queue.total === 0 ? <p className="empty-state">A fila está vazia.</p> : (
            <div className="progress-layout progress-layout--stack">
              <RingGauge share={queue.ok / queue.total} label="Concluídas" caption={`${n(queue.ok)} de ${n(queue.total)} moléculas`} color={statusColor.good} />
              <ProportionBar segments={queueSegments(queue)} total={queue.total} label="Estados da fila" />
            </div>
          ))}
        </section>
        <section className="panel" aria-label="Semi-Imperium">
          <SectionHeading eyebrow="Pipeline canônico" title="Store Semi-Imperium">Registro de execuções e cálculos com proveniência.</SectionHeading>
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

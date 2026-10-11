import type { DatasetProgress } from './api'
import { CountTable, ErrorNotice, Section } from './Resource'
import { MetricCard, MoleculeIllustration, StoryNote } from './StoryVisuals'
import { useApi } from './useApi'

const statusColors: Record<string, string> = {
  OK: '#197e75', Pending: '#a5adaf', Selected: '#6790a0', Running: '#5c77a9',
  Rerun: '#d9a456', Skip: '#c7b1a1', Failed: '#c06459',
}

function StatusBar({ total, status }: { total: number; status: Record<string, number> }) {
  return (
    <div className="status-bar" role="img" aria-label={`Distribuição dos estados do dataset. ${Object.entries(status).map(([name, count]) => `${name}: ${count}`).join('; ')}`}>
      {Object.entries(status).filter(([, count]) => count > 0).map(([name, count]) => (
        <span key={name} title={`${name}: ${count.toLocaleString('pt-BR')}`} style={{ width: `${total ? (count / total) * 100 : 0}%`, backgroundColor: statusColors[name] ?? '#798694' }} />
      ))}
    </div>
  )
}

export function ProgressView() {
  const { state, reload } = useApi<DatasetProgress>('/api/dataset/progress')
  const legacy = state.status === 'ready' ? state.data.legacy : null
  const store = state.status === 'ready' ? state.data.semi_imperium : null
  const complete = legacy?.status.OK ?? 0
  const completion = legacy && legacy.total > 0 ? Math.round((complete / legacy.total) * 100) : null

  return (
    <div className="view">
      <div className="hero hero--progress">
        <div className="hero__copy">
          <span className="eyebrow">01 / Panorama do projeto</span>
          <h1>De moléculas a<br /><em>conhecimento confiável.</em></h1>
          <p>Três perguntas guiam a leitura: o que foi calculado, o que passou pela checagem de qualidade e onde há medidas para comparar.</p>
          <button type="button" className="refresh-button" onClick={reload}>↻ <span>Atualizar dados</span></button>
        </div>
        <div className="hero__art"><MoleculeIllustration /><span>Uma representação visual da rede molecular</span></div>
      </div>

      {state.status === 'loading' && <p role="status" className="loading-state">Carregando progresso…</p>}
      {state.status === 'error' && <ErrorNotice message={state.message} retry={reload} />}
      {legacy && store && (
        <>
          <div className="section-intro"><span className="eyebrow eyebrow--dark">O panorama em números</span><h2>Onde estamos agora</h2><p>{legacy.total > 0 ? `${complete.toLocaleString('pt-BR')} de ${legacy.total.toLocaleString('pt-BR')} registros têm estado OK; ${(legacy.status.Pending ?? 0).toLocaleString('pt-BR')} aguardam processamento (Pending).` : 'Ainda não há registros no CSV PM7.'} O store Semi-Imperium é apresentado separadamente.</p></div>
          <div className="metric-grid">
            <MetricCard label="Registros no dataset" value={legacy.total.toLocaleString('pt-BR')} note="Linhas registradas no CSV PM7" accent />
            <MetricCard label="Cálculo com estado OK" value={complete.toLocaleString('pt-BR')} note={completion === null ? 'Sem linhas para calcular a proporção' : `${completion}% das ${legacy.total.toLocaleString('pt-BR')} linhas`} />
            <MetricCard label="Qualidade CBS marcada OK" value={(legacy.cbs_quality_flag.OK ?? 0).toLocaleString('pt-BR')} note="Flag de qualidade, independente do estado do cálculo" />
          </div>

          <Section title="Dataset PM7">
            <div className="section-heading"><div><span className="eyebrow eyebrow--dark">Etapa 01 / processamento</span><h3>Como as linhas estão distribuídas?</h3><p>O comprimento de cada faixa acompanha a quantidade de linhas em cada estado.</p></div><span className="badge">Protocolo legado: {legacy.protocol}</span></div>
            {legacy.total === 0 ? <p className="empty-state">Nenhuma linha no dataset.</p> : (
              <>
                <StatusBar total={legacy.total} status={legacy.status} />
                <div className="status-legend">{Object.entries(legacy.status).map(([name, count]) => <div key={name}><i style={{ backgroundColor: statusColors[name] ?? '#798694' }} /><span>{name}</span><strong>{count.toLocaleString('pt-BR')}</strong></div>)}</div>
                <div className="detail-grid">
                  <div><h3>Estado do cálculo</h3><p className="subtle">Quantas linhas estão em cada etapa.</p><CountTable counts={legacy.status} label="Estado" /></div>
                  <div><h3>Qualidade da referência CBS</h3><p className="subtle">Uma classificação separada do andamento.</p><CountTable counts={legacy.cbs_quality_flag} label="Flag" /></div>
                </div>
              </>
            )}
          </Section>

          <div className="story-grid">
            <StoryNote number="01" title="O que é PM7?">Um cálculo computacional mais rápido, usado aqui como ponto de partida. Esta tela mostra o protocolo legado crest_pm7.</StoryNote>
            <StoryNote number="02" title="Por que olhar a qualidade CBS?">CBS é uma referência calculada mais rigorosa. Sua flag de qualidade ajuda a identificar quais resultados exigem atenção.</StoryNote>
          </div>

          <Section title="Semi-Imperium">
            <div className="section-heading"><div><span className="eyebrow eyebrow--dark">Etapa 02 / novo registro</span><h3>Um acompanhamento separado</h3><p>Estas contagens vêm do store Semi-Imperium e não entram no total do CSV legado.</p></div></div>
            {!store.exists ? <p className="empty-state">Store ainda não existe.</p> : (
              <>
                <div className="metric-grid metric-grid--compact">
                  <MetricCard label="Moléculas" value={store.molecules.toLocaleString('pt-BR')} note="Identidades registradas" />
                  <MetricCard label="Execuções" value={store.runs.toLocaleString('pt-BR')} note="Rodadas registradas" />
                  <MetricCard label="Cálculos" value={store.calculations.toLocaleString('pt-BR')} note="Trabalhos registrados" />
                </div>
                {store.calculations === 0 ? <p className="empty-state">Nenhum cálculo registrado.</p> : <div className="store-table"><h3>Estado mais recente dos cálculos</h3><CountTable counts={store.states} label="Estado mais recente" /></div>}
              </>
            )}
          </Section>
        </>
      )}
    </div>
  )
}

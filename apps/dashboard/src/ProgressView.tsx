import type { DatasetProgress } from './api'
import { CountTable, ErrorNotice, Section } from './Resource'
import { useApi } from './useApi'

const colors: Record<string, string> = {
  OK: '#2d8659',
  Pending: '#7790a6',
  Selected: '#3d74a5',
  Running: '#7256b0',
  Rerun: '#d3912e',
  Skip: '#bd6464',
  Failed: '#af4545',
}

function StatusBar({ total, status }: { total: number; status: Record<string, number> }) {
  return (
    <div className="status-bar" role="img" aria-label="Distribuição dos estados do dataset">
      {Object.entries(status).map(([name, count]) => (
        <span
          key={name}
          title={`${name}: ${count.toLocaleString('pt-BR')}`}
          style={{ width: `${(count / total) * 100}%`, backgroundColor: colors[name] ?? '#596d7d' }}
        />
      ))}
    </div>
  )
}

export function ProgressView() {
  const { state, reload } = useApi<DatasetProgress>('/api/dataset/progress')

  return (
    <div className="view">
      <div className="view-heading">
        <div>
          <h1>Progresso do dataset</h1>
          <p>Contagens atuais do CSV legado e do store Semi-Imperium.</p>
        </div>
        <button type="button" onClick={reload}>Atualizar dados</button>
      </div>
      {state.status === 'loading' && <p role="status">Carregando progresso…</p>}
      {state.status === 'error' && <ErrorNotice message={state.message} retry={reload} />}
      {state.status === 'ready' && (
        <>
          <Section title="Dataset PM7">
            <p><span className="badge">Protocolo legado: {state.data.legacy.protocol}</span></p>
            <p className="metric">{state.data.legacy.total.toLocaleString('pt-BR')} <small>linhas</small></p>
            {state.data.legacy.total === 0 ? (
              <p>Nenhuma linha no dataset.</p>
            ) : (
              <>
                <StatusBar total={state.data.legacy.total} status={state.data.legacy.status} />
                <div className="columns">
                  <div><h3>Estado do cálculo</h3><CountTable counts={state.data.legacy.status} label="Estado" /></div>
                  <div><h3>Qualidade CBS</h3><CountTable counts={state.data.legacy.cbs_quality_flag} label="Flag" /></div>
                </div>
              </>
            )}
          </Section>
          <Section title="Semi-Imperium">
            {!state.data.semi_imperium.exists ? (
              <p>Store ainda não existe.</p>
            ) : (
              <>
                <div className="metrics">
                  <p><strong>{state.data.semi_imperium.molecules.toLocaleString('pt-BR')}</strong> moléculas</p>
                  <p><strong>{state.data.semi_imperium.runs.toLocaleString('pt-BR')}</strong> execuções</p>
                  <p><strong>{state.data.semi_imperium.calculations.toLocaleString('pt-BR')}</strong> cálculos</p>
                </div>
                {state.data.semi_imperium.calculations === 0 ? (
                  <p>Nenhum cálculo registrado.</p>
                ) : (
                  <CountTable counts={state.data.semi_imperium.states} label="Estado mais recente" />
                )}
              </>
            )}
          </Section>
        </>
      )}
    </div>
  )
}

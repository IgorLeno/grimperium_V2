import { lazy, Suspense, useState } from 'react'
import type { Health } from './api'
import { SummaryView } from './SummaryView'
import type { ViewKey } from './SummaryView'
import { useApi } from './useApi'

const ExperimentalView = lazy(() => import('./ExperimentalView').then((module) => ({ default: module.ExperimentalView })))
const CbsView = lazy(() => import('./CbsView').then((module) => ({ default: module.CbsView })))
const Pm7View = lazy(() => import('./Pm7View').then((module) => ({ default: module.Pm7View })))
const RoutinesView = lazy(() => import('./RoutinesView').then((module) => ({ default: module.RoutinesView })))

const TABS: { key: ViewKey; number: string; label: string; hint: string }[] = [
  { key: 'summary', number: '01', label: 'Resumo', hint: 'Visão geral do projeto' },
  { key: 'experimental', number: '02', label: 'Experimental', hint: 'Medidas de laboratório' },
  { key: 'cbs', number: '03', label: 'CBS', hint: 'Referência de alta precisão' },
  { key: 'pm7', number: '04', label: 'PM7', hint: 'Cálculo rápido (legado)' },
  { key: 'routines', number: '05', label: 'Rotinas', hint: 'Painel de controle' },
]

export function App() {
  const [view, setView] = useState<ViewKey>('summary')
  const health = useApi<Health>('/api/health')
  const navigate = (next: ViewKey) => { setView(next); window.scrollTo?.({ top: 0 }) }
  return (
    <>
      <header className="app-header">
        <div className="brand"><span className="brand-mark" aria-hidden="true">✳</span><span><strong>Grimperium</strong><small>Atlas de dados moleculares</small></span></div>
        <nav aria-label="Telas">
          {TABS.map((tab) => (
            <button key={tab.key} type="button" title={tab.hint} aria-current={view === tab.key ? 'page' : undefined} onClick={() => navigate(tab.key)}>
              <small aria-hidden="true">{tab.number}</small>{tab.label}
            </button>
          ))}
        </nav>
        <div className="api-health" role="status">
          {health.state.status === 'loading' && 'Conectando à API…'}
          {health.state.status === 'ready' && (
            health.state.data.status === 'ok' && health.state.data.read_only
              ? 'API ativa · somente leitura'
              : 'Resposta inesperada da API'
          )}
          {health.state.status === 'error' && (
            <>
              <span title={health.state.message}>API indisponível</span>
              <button type="button" onClick={health.reload}>Reconectar</button>
            </>
          )}
        </div>
      </header>
      <main>
        <Suspense fallback={<p role="status" className="loading-state">Carregando tela…</p>}>
          {view === 'summary' && <SummaryView onNavigate={navigate} />}
          {view === 'experimental' && <ExperimentalView />}
          {view === 'cbs' && <CbsView />}
          {view === 'pm7' && <Pm7View />}
          {view === 'routines' && <RoutinesView />}
        </Suspense>
      </main>
      <footer className="app-footer">Grimperium · dados lidos dos arquivos locais em <code>data/</code> · o painel apenas dispara processos; os cálculos rodam no back-end</footer>
    </>
  )
}

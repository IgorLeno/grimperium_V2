import { lazy, Suspense, useState } from 'react'
import type { Health } from './api'
import { ProgressView } from './ProgressView'
import { useApi } from './useApi'

const ExperimentalView = lazy(() => import('./ExperimentalView').then((module) => ({ default: module.ExperimentalView })))

export function App() {
  const [view, setView] = useState<'progress' | 'experimental'>('progress')
  const health = useApi<Health>('/api/health')
  return (
    <>
      <header className="app-header">
        <strong>Grimperium</strong>
        <nav aria-label="Telas">
          <button type="button" aria-current={view === 'progress' ? 'page' : undefined} onClick={() => setView('progress')}>Progresso</button>
          <button type="button" aria-current={view === 'experimental' ? 'page' : undefined} onClick={() => setView('experimental')}>Experimental</button>
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
      <main>{view === 'progress' ? <ProgressView /> : (
        <Suspense fallback={<p role="status">Carregando tela Experimental…</p>}>
          <ExperimentalView />
        </Suspense>
      )}</main>
    </>
  )
}

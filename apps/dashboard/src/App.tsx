import { Suspense, useState } from 'react'
import type { Health } from './api'
import { DEFAULT_MODULE, firstView, MODULES, topTabOf } from './modules'
import { useApi } from './useApi'

export function App() {
  const [moduleKey, setModuleKey] = useState(DEFAULT_MODULE.key)
  const [view, setView] = useState(firstView(DEFAULT_MODULE.tabs[0]))
  const health = useApi<Health>('/api/health')
  const active = MODULES.find((item) => item.key === moduleKey) ?? DEFAULT_MODULE
  const topTab = topTabOf(active, view)
  const navigate = (next: string) => { setView(next); window.scrollTo?.({ top: 0 }) }
  const selectModule = (key: string) => {
    if (key === active.key) return
    const next = MODULES.find((item) => item.key === key)
    if (!next) return
    setModuleKey(next.key)
    navigate(firstView(next.tabs[0]))
  }
  return (
    <>
      <header className="app-header">
        <div className="brand"><span className="brand-mark" aria-hidden="true">Δ</span><span><strong>Delta Chemistry</strong><small>Métodos rápidos, corrigidos por IA</small></span></div>
        <nav aria-label="Telas">
          {active.tabs.map((tab) => (
            <button key={tab.key} type="button" title={tab.hint} aria-current={topTab.key === tab.key ? 'page' : undefined}
              onClick={() => { if (topTab.key !== tab.key) navigate(firstView(tab)) }}>
              {tab.label}{tab.children && <small aria-hidden="true">▾</small>}
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
      <div className="app-shell">
        <aside className="module-sidebar" aria-label="Módulos">
          <p className="module-sidebar__title">Propriedade</p>
          <nav aria-label="Módulos de propriedade">
            {MODULES.map((item) => (
              <button key={item.key} type="button" aria-current={item.key === active.key ? 'true' : undefined} onClick={() => selectModule(item.key)}>
                <strong>{item.name}</strong>
                <span>{item.property}</span>
                <small>{item.description}</small>
              </button>
            ))}
          </nav>
          <p className="module-sidebar__note">Cada módulo corrige um método rápido para uma propriedade. Novos módulos aparecem aqui.</p>
        </aside>
        <main>
          {topTab.children && (
            <nav className="submenu" aria-label={`Submenu de ${topTab.label}`}>
              <span className="submenu__label">{topTab.label}</span>
              {topTab.children.map((child) => (
                <button key={child.key} type="button" title={child.hint} aria-current={view === child.key ? 'page' : undefined} onClick={() => navigate(child.key)}>{child.label}</button>
              ))}
            </nav>
          )}
          <Suspense fallback={<p role="status" className="loading-state">Carregando tela…</p>}>
            <active.Component view={view} onNavigate={navigate} />
          </Suspense>
        </main>
      </div>
      <footer className="app-footer">Delta Chemistry · dados lidos dos arquivos locais em <code>data/</code> · o painel apenas dispara processos; os cálculos rodam no back-end</footer>
    </>
  )
}

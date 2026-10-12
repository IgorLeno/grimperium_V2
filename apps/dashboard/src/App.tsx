import { Suspense, useState } from 'react'
import type { Health } from './api'
import { DEFAULT_MODULE, firstView, MODULES, topTabOf } from './modules'
import { useApi } from './useApi'

const SIDEBAR_KEY = 'grimperium.sidebar.open'

function readSidebarOpen(): boolean {
  try { return window.localStorage.getItem(SIDEBAR_KEY) !== 'false' } catch { return true }
}

export function App() {
  const [sidebarOpen, setSidebarOpen] = useState(readSidebarOpen)
  const toggleSidebar = () => {
    const next = !sidebarOpen
    setSidebarOpen(next)
    try { window.localStorage.setItem(SIDEBAR_KEY, String(next)) } catch { /* preferência não persistida */ }
  }
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
        <button type="button" className="sidebar-toggle" aria-label={sidebarOpen ? 'Ocultar menu lateral' : 'Mostrar menu lateral'}
          aria-expanded={sidebarOpen} aria-controls="module-sidebar" title={sidebarOpen ? 'Ocultar menu lateral' : 'Mostrar menu lateral'} onClick={toggleSidebar}>
          <svg viewBox="0 0 20 20" width="18" height="18" aria-hidden="true" focusable="false">
            <rect x="2.5" y="3.5" width="15" height="13" rx="3" fill="none" stroke="currentColor" strokeWidth="1.5" />
            <path d="M7.5 3.5v13" stroke="currentColor" strokeWidth="1.5" />
            {sidebarOpen && <rect x="3.25" y="4.25" width="3.5" height="11.5" rx="2" fill="currentColor" opacity=".55" />}
          </svg>
        </button>
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
      <div className={`app-shell${sidebarOpen ? '' : ' app-shell--collapsed'}`}>
        {sidebarOpen && (
          <aside id="module-sidebar" className="module-sidebar" aria-label="Módulos">
            <div className="module-sidebar__head">
              <p className="module-sidebar__title">Módulos</p>
              <span className="module-sidebar__count" aria-hidden="true">{MODULES.length}</span>
            </div>
            <nav aria-label="Módulos de propriedade">
              {MODULES.map((item) => (
                <button key={item.key} type="button" title={item.description} aria-current={item.key === active.key ? 'true' : undefined} onClick={() => selectModule(item.key)}>
                  <span className="module-sidebar__mark" aria-hidden="true">{item.name.charAt(0)}</span>
                  <span className="module-sidebar__text">
                    <strong>{item.name}</strong>
                    <span>{item.property}</span>
                  </span>
                </button>
              ))}
            </nav>
            <p className="module-sidebar__hint">Novos métodos aparecem aqui, um módulo por propriedade.</p>
          </aside>
        )}
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

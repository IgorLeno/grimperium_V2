import { lazy } from 'react'
import { AboutView } from './AboutView'
import { ResultsView } from './ResultsView'
import type { ModuleViewProps } from './modules'

const RoutinesView = lazy(() => import('./RoutinesView').then((module) => ({ default: module.RoutinesView })))

/** Telas do hub Delta Chemistry: Rotinas, Resumos e Sobre. */
export function DeltaModule({ view, onNavigate }: ModuleViewProps) {
  switch (view) {
    case 'results': return <ResultsView onNavigate={onNavigate} />
    case 'about': return <AboutView />
    default: return <RoutinesView eyebrow="Rotinas · Painel de controle" />
  }
}

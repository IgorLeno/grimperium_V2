import { lazy } from 'react'
import { SummaryView } from './SummaryView'
import type { ViewKey } from './SummaryView'
import type { ModuleViewProps } from './modules'

const ExperimentalView = lazy(() => import('./ExperimentalView').then((module) => ({ default: module.ExperimentalView })))
const CbsView = lazy(() => import('./CbsView').then((module) => ({ default: module.CbsView })))
const Pm7View = lazy(() => import('./Pm7View').then((module) => ({ default: module.Pm7View })))
const TrainView = lazy(() => import('./TrainView').then((module) => ({ default: module.TrainView })))
const RoutinesView = lazy(() => import('./RoutinesView').then((module) => ({ default: module.RoutinesView })))

/** Telas do módulo Grimperium (entalpia de formação no estado padrão). */
export function GrimperiumModule({ view, onNavigate }: ModuleViewProps) {
  switch (view as ViewKey) {
    case 'experimental': return <ExperimentalView />
    case 'cbs': return <CbsView />
    case 'pm7': return <Pm7View />
    case 'routines': return <RoutinesView />
    case 'train': return <TrainView />
    default: return <SummaryView onNavigate={onNavigate} />
  }
}

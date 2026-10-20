import { lazy } from 'react'
import { ComingSoonView } from './ComingSoonView'
import { SummaryView } from './SummaryView'
import type { ViewKey } from './SummaryView'
import type { ModuleViewProps } from './modules'

const ExperimentalView = lazy(() => import('./ExperimentalView').then((module) => ({ default: module.ExperimentalView })))
const CbsView = lazy(() => import('./CbsView').then((module) => ({ default: module.CbsView })))
const TrainView = lazy(() => import('./TrainView').then((module) => ({ default: module.TrainView })))
const TrainingView = lazy(() => import('./TrainingView').then((module) => ({ default: module.TrainingView })))

/** Telas do módulo Delta Entalpia: calcular PM7, treinar modelo, prever ΔH e dados. */
export function DeltaEnthalpyModule({ view, onNavigate }: ModuleViewProps) {
  switch (view as ViewKey) {
    case 'experimental': return <ExperimentalView />
    case 'cbs': return <CbsView />
    case 'training': return <TrainingView />
    case 'learning': return <TrainView />
    case 'calc': return (
      <ComingSoonView eyebrow="Prever ΔH · Grimperium" title="Preveja uma molécula" accent="pelo nome."
        lead="Você dará o nome; a rotina reaproveita o valor já calculado ou roda o protocolo e aplica a correção do modelo treinado."
        depends="Depende de Treinar modelo: o modelo corrigido ainda não existe." />
    )
    default: return <SummaryView onNavigate={onNavigate} />
  }
}

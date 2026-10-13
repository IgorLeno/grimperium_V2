import { AboutView } from './AboutView'
import { ResultsView } from './ResultsView'
import type { ModuleViewProps } from './modules'

/** Telas do hub Delta Chemistry: Resumos e Sobre. */
export function DeltaModule({ view }: ModuleViewProps) {
  switch (view) {
    case 'about': return <AboutView />
    default: return <ResultsView />
  }
}

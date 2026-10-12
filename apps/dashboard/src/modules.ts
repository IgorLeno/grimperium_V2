import type { ComponentType } from 'react'
import { GrimperiumModule } from './GrimperiumModule'

/** Contrato que toda tela de módulo recebe do shell do Delta Chemistry. */
export interface ModuleViewProps {
  view: string
  onNavigate: (view: string) => void
}

/** Item do menu superior; `children` forma o submenu (o item pai só agrupa). */
export interface ModuleTab {
  key: string
  label: string
  hint: string
  children?: ModuleTab[]
}

/**
 * Um módulo é um método rápido corrigido por Δ-learning para UMA propriedade.
 * Para adicionar um módulo: crie o componente de telas (como `GrimperiumModule`)
 * e acrescente uma entrada em `MODULES`; a barra lateral e o menu superior
 * passam a reconhecê-lo sem outras alterações no shell.
 */
export interface ModuleDef {
  key: string
  name: string
  property: string
  description: string
  tabs: ModuleTab[]
  Component: ComponentType<ModuleViewProps>
}

export const MODULES: ModuleDef[] = [
  {
    key: 'grimperium',
    name: 'Grimperium',
    property: 'Entalpia de formação (estado padrão)',
    description: 'PM7 corrigido por Δ-learning a partir de cálculos CBS.',
    tabs: [
      { key: 'summary', label: 'Resumo', hint: 'Visão geral do módulo' },
      {
        key: 'data', label: 'Dados', hint: 'Tabelas usadas pelo módulo',
        children: [
          { key: 'experimental', label: 'Experimental', hint: 'Medidas de laboratório' },
          { key: 'cbs', label: 'CBS', hint: 'Referência de alta precisão' },
          { key: 'pm7', label: 'PM7', hint: 'Cálculo rápido (legado)' },
        ],
      },
      { key: 'routines', label: 'Calcular', hint: 'Monta rotinas de cálculo' },
      { key: 'train', label: 'Treinar', hint: 'Treino do modelo de correção Δ' },
    ],
    Component: GrimperiumModule,
  },
]

export const DEFAULT_MODULE = MODULES[0]

/** Primeira tela (folha) de um item de menu. */
export function firstView(tab: ModuleTab): string {
  return tab.children?.[0]?.key ?? tab.key
}

/** Item de topo que contém a tela `view` (ela mesma ou uma filha). */
export function topTabOf(module: ModuleDef, view: string): ModuleTab {
  return module.tabs.find((tab) => tab.key === view || tab.children?.some((child) => child.key === view)) ?? module.tabs[0]
}

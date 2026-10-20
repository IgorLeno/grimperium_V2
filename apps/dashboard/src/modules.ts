import type { ComponentType } from 'react'
import { DeltaModule } from './DeltaModule'
import { DeltaEnthalpyModule } from './DeltaEnthalpyModule'

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
 * Entrada da barra lateral. O hub (Delta Chemistry) e cada módulo de
 * propriedade têm o mesmo formato; o menu superior mostra só as abas do item
 * ativo, então as telas de um módulo não aparecem fora dele.
 */
export interface ModuleDef {
  key: string
  name: string
  /** Texto curto sob o nome na barra lateral. */
  property: string
  description: string
  tabs: ModuleTab[]
  Component: ComponentType<ModuleViewProps>
}

/** Tela principal: resumos e a explicação do projeto. */
export const HUB: ModuleDef = {
  key: 'delta-chemistry',
  name: 'Delta Chemistry',
  property: 'Resumos e sobre',
  description: 'Hub de rotinas de cálculo de propriedades com métodos rápidos corrigidos por machine learning.',
  tabs: [
    { key: 'results', label: 'Resumos', hint: 'Resultados obtidos pelas rotinas' },
    { key: 'about', label: 'Sobre', hint: 'O que é o Delta Chemistry' },
  ],
  Component: DeltaModule,
}

/**
 * Um módulo é um método rápido corrigido por Δ-learning para UMA propriedade.
 * Para adicionar um módulo: crie o componente de telas (como `DeltaEnthalpyModule`)
 * e acrescente uma entrada aqui; a barra lateral e o menu superior passam a
 * reconhecê-lo sem outras alterações no shell.
 */
export const MODULES: ModuleDef[] = [
  {
    key: 'delta-enthalpy',
    name: 'Delta Entalpia',
    property: 'Entalpia de formação (estado padrão)',
    description: 'PM7 corrigido por Δ-learning a partir de cálculos CBS: treinamento, Delta Learning e cálculo.',
    tabs: [
      { key: 'summary', label: 'Resumo', hint: 'Visão geral do módulo' },
      {
        key: 'data', label: 'Dados', hint: 'Tabelas usadas pelo módulo',
        children: [
          { key: 'experimental', label: 'Experimental', hint: 'Medidas de laboratório' },
          { key: 'cbs', label: 'CBS', hint: 'Referência de alta precisão' },
        ],
      },
      { key: 'training', label: 'Treinamento', hint: 'Semi-Imperium: calcula as moléculas PM7 do treino' },
      { key: 'learning', label: 'Delta Learning', hint: 'Treina o modelo de correção Δ' },
      { key: 'calc', label: 'Cálculo', hint: 'Grimperium: calcula uma molécula pelo nome' },
    ],
    Component: DeltaEnthalpyModule,
  },
]

export const DEFAULT_MODULE = HUB

/** Hub + módulos, na ordem em que o shell os procura. */
export const ALL_MODULES: ModuleDef[] = [HUB, ...MODULES]

/** Primeira tela (folha) de um item de menu. */
export function firstView(tab: ModuleTab): string {
  return tab.children?.[0]?.key ?? tab.key
}

/** Item de topo que contém a tela `view` (ela mesma ou uma filha). */
export function topTabOf(module: ModuleDef, view: string): ModuleTab {
  return module.tabs.find((tab) => tab.key === view || tab.children?.some((child) => child.key === view)) ?? module.tabs[0]
}

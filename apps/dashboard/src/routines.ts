import type { RoutineCatalog, RoutineProgram } from './api'

export type ParamValue = string | number | boolean | null
export type ParamSet = Record<string, ParamValue>

/** The contract a future job backend will receive (version bumps on breaking change). */
export interface RoutineSpec {
  version: 1
  id: string
  name: string
  description: string
  property: string
  input: { mode: string; params: ParamSet; molecules: string[] }
  steps: string[]
  programs: Record<string, ParamSet>
  outputs: { kinds: string[]; params: Record<string, ParamSet> }
  created_at: string
  updated_at: string
}

export type Preset = 'full' | 'direct_pm7' | 'custom'

export const PRESETS: Record<Exclude<Preset, 'custom'>, { label: string; steps: string[]; text: string }> = {
  full: { label: 'Completo (Semi-Imperium)', steps: ['rdkit', 'crest', 'confpass', 'mopac'], text: 'Busca conformações e escolhe as melhores antes do PM7. Mais lento, mais confiável.' },
  direct_pm7: { label: 'Direto no PM7', steps: ['rdkit', 'mopac'], text: 'Uma única estrutura 3D do RDKit vai direto ao MOPAC. Rápido, sem busca conformacional.' },
}

export function defaults(program: RoutineProgram): ParamSet {
  return Object.fromEntries(program.parameters.map((parameter) => [parameter.key, parameter.default]))
}

export function newRoutine(catalog: RoutineCatalog): RoutineSpec {
  const now = new Date().toISOString()
  const input = catalog.input_modes[0]
  return {
    version: 1,
    id: `routine-${Date.now().toString(36)}`,
    name: '',
    description: '',
    property: catalog.properties.find((item) => item.available)?.key ?? '',
    input: { mode: input?.key ?? 'smiles', params: input ? defaults(input) : {}, molecules: [] },
    steps: [...PRESETS.full.steps],
    programs: Object.fromEntries(catalog.programs.map((program) => [program.key, defaults(program)])),
    outputs: { kinds: ['spreadsheet'], params: Object.fromEntries(catalog.outputs.map((output) => [output.key, defaults(output)])) },
    created_at: now,
    updated_at: now,
  }
}

export function presetOf(steps: string[]): Preset {
  for (const [key, preset] of Object.entries(PRESETS)) {
    if (preset.steps.join() === steps.join()) return key as Preset
  }
  return 'custom'
}

/** Ordering rules of the Semi-Imperium pipeline, in plain language. */
export function pipelineProblems(steps: string[]): string[] {
  const problems: string[] = []
  const at = (key: string) => steps.indexOf(key)
  if (!steps.includes('rdkit')) problems.push('O RDKit é obrigatório: ele valida a molécula e gera a estrutura 3D inicial.')
  if (!steps.includes('mopac')) problems.push('O MOPAC é obrigatório para calcular a entalpia de formação.')
  if (steps.includes('confpass') && !steps.includes('crest')) problems.push('A seleção de conformações precisa do CREST antes dela.')
  if (at('rdkit') > 0) problems.push('O RDKit deve ser o primeiro passo.')
  if (steps.includes('mopac') && at('mopac') !== steps.length - 1) problems.push('O MOPAC deve ser o último passo.')
  if (steps.includes('crest') && steps.includes('confpass') && at('confpass') < at('crest')) problems.push('A seleção de conformações vem depois do CREST.')
  return problems
}

export function routineProblems(spec: RoutineSpec): { step: number; text: string }[] {
  const problems: { step: number; text: string }[] = []
  if (!spec.name.trim()) problems.push({ step: 1, text: 'Dê um nome à rotina.' })
  if (!spec.property) problems.push({ step: 2, text: 'Escolha a propriedade a calcular.' })
  if (spec.input.molecules.length === 0) problems.push({ step: 3, text: 'Informe ao menos uma molécula.' })
  for (const text of pipelineProblems(spec.steps)) problems.push({ step: 3, text })
  if (spec.outputs.kinds.length === 0) problems.push({ step: 5, text: 'Escolha ao menos um tipo de saída.' })
  return problems
}

export function parseMolecules(text: string): string[] {
  return text.split('\n').map((line) => line.trim()).filter(Boolean)
}

const STORAGE_KEY = 'grimperium.routines.v1'

/** Routines live in this browser until the job backend stores them. */
export function loadRoutines(): RoutineSpec[] {
  try {
    const raw = window.localStorage.getItem(STORAGE_KEY)
    const parsed: unknown = raw ? JSON.parse(raw) : []
    return Array.isArray(parsed) ? (parsed as RoutineSpec[]).filter((item) => item && item.version === 1) : []
  } catch {
    return []
  }
}

export function saveRoutines(routines: RoutineSpec[]): boolean {
  try {
    window.localStorage.setItem(STORAGE_KEY, JSON.stringify(routines))
    return true
  } catch {
    return false
  }
}

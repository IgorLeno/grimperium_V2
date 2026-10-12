export interface Health {
  status: string
  read_only: boolean
}

export interface DatasetProgress {
  legacy: {
    protocol: string
    path: string
    total: number
    status: Record<string, number>
    cbs_quality_flag: Record<string, number>
  }
  semi_imperium: {
    path: string
    exists: boolean
    molecules: number
    runs: number
    calculations: number
    states: Record<string, number>
  }
}

export interface ExperimentalSummary {
  path: string
  rows: number
  by_source: Record<string, number>
  by_source_family: Record<string, number>
  flags: Record<string, number>
  uncertainty_kind: Record<string, number>
}

export interface ComparisonPoint {
  exp_id: string
  mol_id: string
  smiles: string
  nheavy: number
  source: string
  H298_exp: number
  uncertainty: number | null
  H298_cbs: number
  cbs_minus_exp: number
  conflict_flag: boolean
  low_priority_only: boolean
  validation_only: boolean
}

export interface ExperimentalComparison {
  path: string
  units: string
  stats: {
    n: number
    mae: number | null
    bias: number | null
    rmse: number | null
  }
  points: ComparisonPoint[]
}

export async function getJson<T>(path: string, signal: AbortSignal): Promise<T> {
  let response: Response
  try {
    response = await fetch(path, { signal })
  } catch (error) {
    if (signal.aborted) throw error
    throw new Error(`Não foi possível conectar à API (${path}). Verifique se o servidor está ativo.`, { cause: error })
  }

  if (!response.ok) {
    const body: unknown = await response.json().catch(() => null)
    const detail =
      typeof body === 'object' && body !== null && 'detail' in body
        ? String(body.detail)
        : response.statusText || 'erro desconhecido'
    const file =
      typeof body === 'object' && body !== null && 'path' in body
        ? ` (${String(body.path)})`
        : ''
    const hint = response.status === 500 || response.status === 502
      ? ' Verifique se a API está ativa na porta 8000.'
      : ''
    throw new Error(`${path}: ${response.status} — ${detail}${file}${hint}`)
  }
  return (await response.json()) as T
}

export type TableKey = 'exp' | 'cbs' | 'pm7'

export interface HistogramBin {
  start: number
  end: number
  count: number
}

export interface Histogram {
  column: string
  label: string
  unit: string
  n: number
  below: number
  above: number
  mean: number | null
  median: number | null
  p05: number | null
  p95: number | null
  minimum: number | null
  maximum: number | null
  bins: HistogramBin[]
}

export interface SizeStats {
  nheavy: number
  n: number
  mean: number
  mean_abs: number
}

export interface ProfilePair {
  id: string
  x: number
  y: number
  nheavy: number | null
}

export interface DatasetProfile {
  key: TableKey
  title: string
  path: string
  rows: number
  unique_smiles: number | null
  nheavy: Record<string, number>
  composition: Record<'CH' | 'CHO' | 'CHN' | 'CHNO', number>
  categories: Record<string, Record<string, number>>
  histograms: Histogram[]
  by_nheavy: SizeStats[]
  pairs: ProfilePair[]
  pairs_total: number
  missing_columns: string[]
  manifest: Record<string, unknown> | null
}

export type Cell = string | number | boolean | null

export interface TablePage {
  key: TableKey
  columns: string[]
  total: number
  matched: number
  offset: number
  limit: number
  rows: Record<string, Cell>[]
}

export interface RoutineParameter {
  key: string
  label: string
  kind: 'int' | 'float' | 'bool' | 'choice' | 'text'
  default: string | number | boolean | null
  choices: string[]
  unit: string
  help: string
  minimum: number | null
  maximum: number | null
}

export interface RoutineProgram {
  key: string
  label: string
  role: string
  description: string
  executable: string | null
  available: boolean | null
  experimental: boolean
  parameters: RoutineParameter[]
}

export interface RoutineProperty {
  key: string
  label: string
  unit: string
  description: string
  available: boolean
}

export interface RoutineCatalog {
  dispatch_enabled: boolean
  properties: RoutineProperty[]
  input_modes: RoutineProgram[]
  programs: RoutineProgram[]
  outputs: RoutineProgram[]
}

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

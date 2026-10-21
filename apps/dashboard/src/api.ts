export interface Health {
  status: string
  can_dispatch: boolean
}

export interface DatasetProgress {
  queue: TrainingProgress
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

async function requestJson<T>(path: string, init: RequestInit): Promise<T> {
  let response: Response
  try {
    response = await fetch(path, init)
  } catch (error) {
    if (init.signal?.aborted) throw error
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

export function getJson<T>(path: string, signal: AbortSignal): Promise<T> {
  return requestJson<T>(path, { signal })
}

export function postJson<T>(path: string, body: unknown = {}): Promise<T> {
  return requestJson<T>(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  })
}

export type TableKey = 'exp' | 'cbs'

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

/** Moléculas da fila de treinamento por etapa; `total` é o denominador. */
export interface TrainingProgress {
  total: number
  ok: number
  pending: number
  rerun: number
  running: number
  skipped: number
}

export interface TrainingJob {
  state: 'idle' | 'running' | 'stopping' | 'finished' | 'failed' | 'stale'
  pid: number | null
  started_at: string | null
  heartbeat_at: string | null
  finished_at: string | null
  workers: number
  limit: number | null
  done: number
  failed: number
  message: string
  stop_requested: boolean
}

export interface TrainingFailure { mol_id: string; reason: string; at: string }

/** Molécula marcada como em cálculo; `started_at` é quando foi despachada. */
export interface TrainingInFlight { mol_id: string; smiles: string; nheavy: number; started_at: string | null }

export interface TrainingCompleted {
  mol_id: string
  smiles: string
  h298_pm7: number
  delta_kcal_mol: number | null
  elapsed_s: number | null
  conformers_evaluated: number | null
  finished_at: string
}

export interface TrainingActivity {
  in_flight: TrainingInFlight[]
  recent: TrainingCompleted[]
  /** Tempo médio por molécula aceita desde o início da rodada. */
  mean_elapsed_s: number | null
}

export interface TrainingStatus {
  progress: TrainingProgress
  job: TrainingJob
  geometries: number
  failures: TrainingFailure[]
  activity: TrainingActivity
  cpu_count: number
}

/** Uma célula (computador remoto) como a central a vê; nunca o token. */
export interface CellSummary {
  cell_id: string
  name: string
  /** online · late · silent · never_seen · version_mismatch · blocked · revoked */
  state: string
  admission: string
  host: string | null
  cores: number | null
  workers: number | null
  crest_threads: number | null
  program_versions: Record<string, string>
  refusal: string | null
  created_at: string
  registered_at: string | null
  last_seen_at: string | null
  revoked_at: string | null
  held: number
}

export interface CellCreated { cell: CellSummary; token: string }

/**
 * Pede o pacote da célula (zip com token novo) e dispara o download.
 * O zip só existe na memória do navegador até ser salvo.
 */
export async function downloadCellPackage(cellId: string, centralUrl: string): Promise<string> {
  const path = `/api/cells/${encodeURIComponent(cellId)}/package`
  let response: Response
  try {
    response = await fetch(path, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ central_url: centralUrl }),
    })
  } catch (error) {
    throw new Error(`Não foi possível conectar à API (${path}). Verifique se o servidor está ativo.`, { cause: error })
  }
  if (!response.ok) {
    const body: unknown = await response.json().catch(() => null)
    const detail = typeof body === 'object' && body !== null && 'detail' in body
      ? (typeof body.detail === 'string' ? body.detail : 'endereço da central inválido (use http://IP:porta)')
      : response.statusText || 'erro desconhecido'
    throw new Error(`${response.status} — ${detail}`)
  }
  const match = /filename="([^"]+)"/.exec(response.headers.get('Content-Disposition') ?? '')
  const filename = match?.[1] ?? 'celula.zip'
  const url = URL.createObjectURL(await response.blob())
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  document.body.append(link)
  link.click()
  link.remove()
  // Revogar na hora pode cancelar o download em alguns navegadores.
  window.setTimeout(() => URL.revokeObjectURL(url), 60_000)
  return filename
}

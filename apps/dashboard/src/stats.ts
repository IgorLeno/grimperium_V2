import type { ComparisonPoint } from './api'

/** 1 kcal/mol: the conventional "chemical accuracy" target for enthalpies. */
export const CHEMICAL_ACCURACY = 1

export interface ErrorSummary {
  n: number
  mae: number | null
  bias: number | null
  rmse: number | null
  /** Share (0–1) of pairs with |CBS − exp| ≤ CHEMICAL_ACCURACY. */
  withinChemical: number | null
  /** Median of |CBS − exp|. */
  medianAbs: number | null
}

export function median(values: number[]): number | null {
  if (values.length === 0) return null
  const sorted = [...values].sort((a, b) => a - b)
  const middle = Math.floor(sorted.length / 2)
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2
}

export function summarize(points: ComparisonPoint[]): ErrorSummary {
  const n = points.length
  if (n === 0) return { n, mae: null, bias: null, rmse: null, withinChemical: null, medianAbs: null }
  const errors = points.map((point) => point.cbs_minus_exp)
  const abs = errors.map(Math.abs)
  return {
    n,
    mae: abs.reduce((sum, value) => sum + value, 0) / n,
    bias: errors.reduce((sum, value) => sum + value, 0) / n,
    rmse: Math.sqrt(errors.reduce((sum, value) => sum + value ** 2, 0) / n),
    withinChemical: abs.filter((value) => value <= CHEMICAL_ACCURACY).length / n,
    medianAbs: median(abs),
  }
}

export type TierKey = 'excellent' | 'close' | 'attention' | 'outlier'

export interface Tier {
  key: TierKey
  label: string
  range: string
  description: string
  /** Inclusive upper bound of |error| in kcal/mol. */
  max: number
}

export const TIERS: Tier[] = [
  { key: 'excellent', label: 'Excelente', range: '≤ 1', max: CHEMICAL_ACCURACY, description: 'Dentro da precisão química' },
  { key: 'close', label: 'Próximo', range: '1 – 5', max: 5, description: 'Diferença pequena' },
  { key: 'attention', label: 'Atenção', range: '5 – 10', max: 10, description: 'Diferença perceptível' },
  { key: 'outlier', label: 'Discrepante', range: '> 10', max: Infinity, description: 'Merece investigação' },
]

export function tierOf(error: number): Tier {
  const abs = Math.abs(error)
  return TIERS.find((tier) => abs <= tier.max) ?? TIERS[TIERS.length - 1]
}

export function tierCounts(points: ComparisonPoint[]): { tier: Tier; count: number }[] {
  const counts = new Map<TierKey, number>(TIERS.map((tier) => [tier.key, 0]))
  for (const point of points) {
    const key = tierOf(point.cbs_minus_exp).key
    counts.set(key, (counts.get(key) ?? 0) + 1)
  }
  return TIERS.map((tier) => ({ tier, count: counts.get(tier.key) ?? 0 }))
}

export interface HistogramBin {
  start: number
  end: number
  /** Bin centre, used as the x position. */
  mid: number
  label: string
  count: number
}

export function quantile(sorted: number[], q: number): number {
  if (sorted.length === 0) return NaN
  const position = (sorted.length - 1) * q
  const low = Math.floor(position)
  const high = Math.ceil(position)
  return sorted[low] + (sorted[high] - sorted[low]) * (position - low)
}

export interface ClippedHistogram {
  bins: HistogramBin[]
  /** Values left out of the axis range, below and above it. */
  below: number
  above: number
}

/**
 * Fixed-width histogram of signed errors over the central ``1 − 2·clip`` share
 * of the values; bins are aligned to multiples of ``width`` and the tails are
 * counted, not drawn, so a few outliers cannot flatten the chart.
 */
export function histogram(values: number[], width = 2, clip = 0.005): ClippedHistogram {
  if (values.length === 0) return { bins: [], below: 0, above: 0 }
  const sorted = [...values].sort((a, b) => a - b)
  const lo = Math.floor(quantile(sorted, clip) / width) * width
  let hi = Math.ceil(quantile(sorted, 1 - clip) / width) * width
  if (hi <= lo) hi = lo + width
  const bins: HistogramBin[] = []
  for (let start = lo; start < hi; start += width) {
    const end = start + width
    bins.push({ start, end, mid: start + width / 2, label: `${start} a ${end}`, count: 0 })
  }
  let below = 0
  let above = 0
  for (const value of values) {
    if (value < lo) below += 1
    else if (value > hi) above += 1
    else bins[Math.min(bins.length - 1, Math.floor((value - lo) / width))].count += 1
  }
  return { bins, below, above }
}

export interface CumulativePoint {
  threshold: number
  share: number
}

/** Share (0–100) of pairs whose |error| is at most each threshold. */
export function cumulativeAccuracy(points: ComparisonPoint[], step = 0.5): CumulativePoint[] {
  if (points.length === 0) return []
  const abs = points.map((point) => Math.abs(point.cbs_minus_exp)).sort((a, b) => a - b)
  const maximum = Math.ceil(abs[abs.length - 1] / step) * step
  const curve: CumulativePoint[] = []
  let index = 0
  for (let threshold = 0; threshold <= maximum + 1e-9; threshold += step) {
    while (index < abs.length && abs[index] <= threshold + 1e-9) index += 1
    curve.push({ threshold: Number(threshold.toFixed(2)), share: (index / abs.length) * 100 })
  }
  return curve
}

/** Smallest threshold that covers at least ``share`` percent of pairs. */
export function thresholdFor(points: ComparisonPoint[], share: number): number | null {
  if (points.length === 0) return null
  const abs = points.map((point) => Math.abs(point.cbs_minus_exp)).sort((a, b) => a - b)
  const index = Math.min(abs.length - 1, Math.max(0, Math.ceil((share / 100) * abs.length) - 1))
  return abs[index]
}

export interface GroupStats {
  key: string
  n: number
  mae: number
  bias: number
}

function group(points: ComparisonPoint[], keyOf: (point: ComparisonPoint) => string): GroupStats[] {
  const groups = new Map<string, number[]>()
  for (const point of points) {
    const key = keyOf(point)
    const errors = groups.get(key)
    if (errors) errors.push(point.cbs_minus_exp)
    else groups.set(key, [point.cbs_minus_exp])
  }
  return [...groups.entries()].map(([key, errors]) => ({
    key,
    n: errors.length,
    mae: errors.reduce((sum, value) => sum + Math.abs(value), 0) / errors.length,
    bias: errors.reduce((sum, value) => sum + value, 0) / errors.length,
  }))
}

/** Error per source, keeping sources with at least ``minimum`` pairs, lowest MAE first. */
export function bySource(points: ComparisonPoint[], minimum = 5): GroupStats[] {
  return group(points, (point) => point.source)
    .filter((item) => item.n >= minimum)
    .sort((a, b) => a.mae - b.mae || a.key.localeCompare(b.key))
}

export interface SizeStats extends GroupStats {
  nheavy: number
}

/** Error per number of heavy atoms, ascending size. */
export function byHeavyAtoms(points: ComparisonPoint[]): SizeStats[] {
  return group(points, (point) => String(point.nheavy))
    .map((item) => ({ ...item, nheavy: Number(item.key) }))
    .sort((a, b) => a.nheavy - b.nheavy)
}

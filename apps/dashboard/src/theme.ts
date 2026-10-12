/**
 * Chart colors as CSS custom properties, so light and dark mode swap in
 * style.css only. The categorical order was validated for color-vision
 * deficiency (adjacent pairs; the first three also all-pairs) on both
 * surfaces — keep the order when adding series.
 */
export const series = ['var(--series-1)', 'var(--series-2)', 'var(--series-3)', 'var(--series-4)'] as const

/** Meaning-bearing colors; always shown next to an icon or a label. */
export const statusColor = {
  good: 'var(--status-good)',
  warning: 'var(--status-warning)',
  serious: 'var(--status-serious)',
  critical: 'var(--status-critical)',
  neutral: 'var(--status-neutral)',
  info: 'var(--status-info)',
} as const

export const chrome = {
  grid: 'var(--chart-grid)',
  axis: 'var(--chart-axis)',
  muted: 'var(--text-muted)',
  ink: 'var(--text-secondary)',
  surface: 'var(--surface-1)',
} as const

export const axisTick = { fill: 'var(--text-muted)', fontSize: 12 }

export function formatNumber(value: number | null | undefined, digits = 2): string {
  return value === null || value === undefined || Number.isNaN(value)
    ? '—'
    : value.toLocaleString('pt-BR', { maximumFractionDigits: digits })
}

export function formatSigned(value: number | null | undefined, digits = 2): string {
  if (value === null || value === undefined || Number.isNaN(value)) return '—'
  const text = Math.abs(value).toLocaleString('pt-BR', { maximumFractionDigits: digits })
  return value > 0 ? `+${text}` : value < 0 ? `−${text}` : text
}

export function formatPercent(share: number | null | undefined, digits = 0): string {
  return share === null || share === undefined || Number.isNaN(share)
    ? '—'
    : `${(share * 100).toLocaleString('pt-BR', { maximumFractionDigits: digits })}%`
}

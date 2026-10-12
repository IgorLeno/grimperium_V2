import type { ReactNode } from 'react'

interface TooltipProps<T> {
  active?: boolean
  payload?: ReadonlyArray<{ payload?: T }>
}

/** Builds a Recharts ``content`` renderer that shows one row of the hovered datum. */
export function tooltipFor<T>(render: (datum: T) => ReactNode) {
  return function ChartTooltip({ active, payload }: TooltipProps<T>) {
    const datum = payload?.[0]?.payload
    if (!active || datum === undefined) return null
    return <div className="chart-tooltip">{render(datum)}</div>
  }
}

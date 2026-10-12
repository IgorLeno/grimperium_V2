import {
  Area,
  AreaChart,
  Bar,
  BarChart,
  CartesianGrid,
  Cell,
  ComposedChart,
  LabelList,
  Line,
  Pie,
  PieChart,
  ReferenceLine,
  ResponsiveContainer,
  Scatter,
  ScatterChart,
  Tooltip,
  XAxis,
  YAxis,
} from 'recharts'
import type { Histogram } from './api'
import type { ReactNode } from 'react'
import { tooltipFor } from './ChartTooltip'
import { niceTicks } from './helpers'
import type { Segment } from './StoryVisuals'
import { axisTick, chrome, formatNumber, formatSigned, series } from './theme'

function TooltipRow({ label, value, color }: { label: string; value: ReactNode; color?: string }) {
  return (
    <div className="chart-tooltip__row">
      <span>{color && <i style={{ backgroundColor: color }} aria-hidden="true" />}{label}</span>
      <strong>{value}</strong>
    </div>
  )
}

const margin = { top: 22, right: 16, bottom: 28, left: 8 }


export interface RefMark {
  value: number
  label: string
  color?: string
}

interface BinDatum { mid: number; start: number; end: number; count: number }

/** Histogram from the API (or from stats.ts), with labelled reference lines. */
export function HistogramChart({ histogram, color = series[0], xLabel, refs = [], height = 260 }: {
  histogram: Pick<Histogram, 'bins' | 'unit'>
  color?: string
  xLabel: string
  refs?: RefMark[]
  height?: number
}) {
  const data: BinDatum[] = histogram.bins.map((bin) => ({ ...bin, mid: (bin.start + bin.end) / 2 }))
  const width = data.length ? data[0].end - data[0].start : 1
  const Content = tooltipFor<BinDatum>((bin) => (
    <>
      <p className="chart-tooltip__title">{formatNumber(bin.start)} a {formatNumber(bin.end)} {histogram.unit}</p>
      <TooltipRow label="Moléculas" value={bin.count.toLocaleString('pt-BR')} color={color} />
    </>
  ))
  return (
    <div className="chart" role="img" aria-label={`Histograma: ${xLabel}`}>
      <ResponsiveContainer width="100%" height={height}>
        <BarChart data={data} margin={margin} barCategoryGap={1}>
          <CartesianGrid vertical={false} stroke={chrome.grid} />
          <XAxis
            dataKey="mid" type="number" domain={[data[0]?.start ?? 0, data[data.length - 1]?.end ?? 1]}
            ticks={niceTicks(data[0]?.start ?? 0, data[data.length - 1]?.end ?? 1)}
            tick={axisTick} tickLine={false} axisLine={{ stroke: chrome.axis }}
            tickFormatter={(value: number) => formatNumber(value, 1)}
            label={{ value: xLabel, position: 'insideBottom', offset: -16, fill: chrome.muted, fontSize: 12 }}
          />
          <YAxis tick={axisTick} tickLine={false} axisLine={false} allowDecimals={false} width={44} />
          <Tooltip content={Content} cursor={{ fill: 'var(--hover-wash)' }} />
          <Bar dataKey="count" fill={color} radius={[4, 4, 0, 0]} maxBarSize={48} />
          {[...refs].sort((a, b) => a.value - b.value).map((mark, index, sorted) => (
            <ReferenceLine
              key={mark.label} x={mark.value} stroke={mark.color ?? 'var(--text-primary)'} strokeWidth={1.5}
              strokeDasharray={mark.value === 0 ? undefined : '4 3'}
              label={{
                value: mark.label, position: 'top', fill: 'var(--text-secondary)', fontSize: 11,
                // Two nearby lines: push the labels apart instead of stacking them.
                dx: sorted.length > 1 ? (index === 0 ? -16 : 16) : 0,
              }}
              ifOverflow="extendDomain"
            />
          ))}
        </BarChart>
      </ResponsiveContainer>
      <span className="sr-only">Largura de cada barra: {formatNumber(width)} {histogram.unit}</span>
    </div>
  )
}

interface CountDatum { name: string; count: number }

/** Vertical columns for an ordered discrete variable (e.g. heavy atoms). */
export function ColumnChart({ counts, xLabel, color = series[0], height = 240, order = 'numeric' }: {
  counts: Record<string, number>
  xLabel: string
  color?: string
  height?: number
  order?: 'numeric' | 'count'
}) {
  const data: CountDatum[] = Object.entries(counts)
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => order === 'numeric' ? Number(a.name) - Number(b.name) : b.count - a.count)
  const Content = tooltipFor<CountDatum>((datum) => (
    <>
      <p className="chart-tooltip__title">{xLabel}: {datum.name}</p>
      <TooltipRow label="Moléculas" value={datum.count.toLocaleString('pt-BR')} color={color} />
    </>
  ))
  return (
    <div className="chart" role="img" aria-label={`${xLabel}. ${data.map((d) => `${d.name}: ${d.count}`).join('; ')}`}>
      <ResponsiveContainer width="100%" height={height}>
        <BarChart data={data} margin={margin}>
          <CartesianGrid vertical={false} stroke={chrome.grid} />
          <XAxis dataKey="name" tick={axisTick} tickLine={false} axisLine={{ stroke: chrome.axis }} interval="preserveStartEnd"
            label={{ value: xLabel, position: 'insideBottom', offset: -16, fill: chrome.muted, fontSize: 12 }} />
          <YAxis tick={axisTick} tickLine={false} axisLine={false} allowDecimals={false} width={44} />
          <Tooltip content={Content} cursor={{ fill: 'var(--hover-wash)' }} />
          <Bar dataKey="count" fill={color} radius={[4, 4, 0, 0]} maxBarSize={36} />
        </BarChart>
      </ResponsiveContainer>
    </div>
  )
}

/** Part-to-whole for ≤ 6 segments, with the headline number in the hole. */
export function DonutChart({ segments, centerValue, centerCaption, label, height = 220 }: {
  segments: Segment[]
  centerValue: string
  centerCaption: string
  label: string
  height?: number
}) {
  const total = segments.reduce((sum, segment) => sum + segment.count, 0)
  const data = segments.filter((segment) => segment.count > 0)
  const Content = tooltipFor<Segment>((segment) => (
    <>
      <p className="chart-tooltip__title">{segment.label}</p>
      <TooltipRow label="Moléculas" value={segment.count.toLocaleString('pt-BR')} color={segment.color} />
      <TooltipRow label="Participação" value={total ? `${Math.round((segment.count / total) * 100)}%` : '—'} />
    </>
  ))
  return (
    <div className="donut">
      <div className="donut__chart" role="img" aria-label={`${label}. ${segments.map((s) => `${s.label}: ${s.count}`).join('; ')}`}>
        <ResponsiveContainer width="100%" height={height}>
          <PieChart>
            <Tooltip content={Content} />
            <Pie data={data} dataKey="count" nameKey="label" innerRadius="64%" outerRadius="92%" paddingAngle={data.length > 1 ? 1.5 : 0}
              stroke="var(--surface-1)" strokeWidth={2} cornerRadius={4} isAnimationActive={false}>
              {data.map((segment) => <Cell key={segment.key} fill={segment.color} />)}
            </Pie>
          </PieChart>
        </ResponsiveContainer>
        <div className="donut__center"><strong>{centerValue}</strong><span>{centerCaption}</span></div>
      </div>
      <ul className="proportion__legend proportion__legend--stack">
        {segments.map((segment) => (
          <li key={segment.key}>
            <i style={{ backgroundColor: segment.color }} aria-hidden="true" />
            <span className="proportion__name">{segment.label}</span>
            <strong>{segment.count.toLocaleString('pt-BR')}</strong>
            <span className="proportion__share">{total ? `${Math.round((segment.count / total) * 100)}%` : '—'}</span>
            {segment.description && <small>{segment.description}</small>}
          </li>
        ))}
      </ul>
    </div>
  )
}

export interface ParityDatum {
  x: number
  y: number
  id: string
  group: number
  detail?: string
  marker?: 'circle' | 'triangle'
}

/** Calculated vs reference, with the y = x line and a ±band of agreement. */
export function ParityChart({ data, groups, xName, yName, band = 5, unit = 'kcal/mol', height = 340 }: {
  data: ParityDatum[]
  groups: { label: string; color: string; marker?: 'circle' | 'triangle'; match: (datum: ParityDatum) => boolean }[]
  xName: string
  yName: string
  band?: number
  unit?: string
  height?: number
}) {
  const values = data.flatMap((datum) => [datum.x, datum.y])
  const lower = values.length ? Math.floor(Math.min(...values) / 50) * 50 : -1
  const upper = values.length ? Math.ceil(Math.max(...values) / 50) * 50 : 1
  const domain: [number, number] = lower === upper ? [lower - 1, upper + 1] : [lower, upper]
  const Content = tooltipFor<ParityDatum>((datum) => (
    <>
      <p className="chart-tooltip__title">{datum.id}</p>
      {datum.detail && <p className="chart-tooltip__detail">{datum.detail}</p>}
      <TooltipRow label={xName} value={`${formatNumber(datum.x)} ${unit}`} />
      <TooltipRow label={yName} value={`${formatNumber(datum.y)} ${unit}`} />
      <TooltipRow label="Diferença" value={`${formatSigned(datum.y - datum.x)} ${unit}`} />
    </>
  ))
  return (
    <div className="chart" role="img" aria-label={`Gráfico de paridade: ${yName} versus ${xName}`}>
      <ResponsiveContainer width="100%" height={height}>
        <ScatterChart margin={{ top: 16, right: 16, bottom: 32, left: 12 }}>
          <CartesianGrid stroke={chrome.grid} />
          <XAxis type="number" dataKey="x" name={xName} domain={domain} ticks={niceTicks(domain[0], domain[1], 5)} tick={axisTick} tickLine={false} axisLine={{ stroke: chrome.axis }}
            label={{ value: `${xName} (${unit})`, position: 'insideBottom', offset: -18, fill: chrome.muted, fontSize: 12 }} />
          <YAxis type="number" dataKey="y" name={yName} domain={domain} ticks={niceTicks(domain[0], domain[1], 5)} tick={axisTick} tickLine={false} axisLine={false} width={52}
            label={{ value: `${yName} (${unit})`, angle: -90, position: 'insideLeft', fill: chrome.muted, fontSize: 12, dy: 50 }} />
          <Tooltip content={Content} cursor={{ stroke: chrome.axis }} />
          <ReferenceLine segment={[{ x: domain[0], y: domain[0] + band }, { x: domain[1] - band, y: domain[1] }]} stroke={chrome.axis} strokeDasharray="4 4" />
          <ReferenceLine segment={[{ x: domain[0] + band, y: domain[0] }, { x: domain[1], y: domain[1] - band }]} stroke={chrome.axis} strokeDasharray="4 4" />
          <ReferenceLine segment={[{ x: domain[0], y: domain[0] }, { x: domain[1], y: domain[1] }]} stroke="var(--text-secondary)" strokeWidth={1.5} />
          {groups.map((group) => {
            const subset = data.filter(group.match)
            return subset.length ? (
              <Scatter key={group.label} name={group.label} data={subset} fill={group.color} fillOpacity={0.72}
                shape={group.marker ?? 'circle'} isAnimationActive={false} />
            ) : null
          })}
        </ScatterChart>
      </ResponsiveContainer>
    </div>
  )
}

export function LegendRow({ items }: { items: { label: string; color?: string; marker?: 'circle' | 'triangle' | 'line' | 'dash' }[] }) {
  return (
    <ul className="legend-row">
      {items.map((item) => (
        <li key={item.label}>
          <span className={`legend-mark legend-mark--${item.marker ?? 'circle'}`} style={item.color ? { color: item.color } : undefined} aria-hidden="true" />
          {item.label}
        </li>
      ))}
    </ul>
  )
}

interface CurveDatum { threshold: number; share: number }

/** "What share of molecules lands within x kcal/mol?" */
export function CumulativeChart({ data, marks, height = 260 }: {
  data: CurveDatum[]
  marks: number[]
  height?: number
}) {
  const shareAt = (x: number) => [...data].reverse().find((point) => point.threshold <= x + 1e-9)?.share ?? 0
  const Content = tooltipFor<CurveDatum>((datum) => (
    <>
      <p className="chart-tooltip__title">Erro de até {formatNumber(datum.threshold, 1)} kcal/mol</p>
      <TooltipRow label="Pares cobertos" value={`${formatNumber(datum.share, 1)}%`} color={series[0]} />
    </>
  ))
  return (
    <div className="chart" role="img" aria-label={`Curva de acerto acumulado. ${marks.map((m) => `até ${m} kcal/mol: ${Math.round(shareAt(m))}%`).join('; ')}`}>
      <ResponsiveContainer width="100%" height={height}>
        <AreaChart data={data} margin={margin}>
          <defs>
            <linearGradient id="cumulative-fill" x1="0" y1="0" x2="0" y2="1">
              <stop offset="0" stopColor={series[0]} stopOpacity={0.35} />
              <stop offset="1" stopColor={series[0]} stopOpacity={0.02} />
            </linearGradient>
          </defs>
          <CartesianGrid vertical={false} stroke={chrome.grid} />
          <XAxis dataKey="threshold" type="number" domain={[0, 'dataMax']} ticks={niceTicks(0, data[data.length - 1]?.threshold ?? 1)} tick={axisTick} tickLine={false} axisLine={{ stroke: chrome.axis }}
            label={{ value: 'Erro absoluto máximo (kcal/mol)', position: 'insideBottom', offset: -16, fill: chrome.muted, fontSize: 12 }} />
          <YAxis domain={[0, 100]} tick={axisTick} tickLine={false} axisLine={false} width={44} tickFormatter={(value: number) => `${value}%`} />
          <Tooltip content={Content} cursor={{ stroke: chrome.axis }} />
          <Area type="stepAfter" dataKey="share" stroke={series[0]} strokeWidth={2} fill="url(#cumulative-fill)" isAnimationActive={false} />
          {marks.map((mark) => <ReferenceLine key={mark} x={mark} stroke={chrome.axis} strokeDasharray="4 3" />)}
        </AreaChart>
      </ResponsiveContainer>
      <ul className="curve-marks">
        {marks.map((mark) => (
          <li key={mark}><strong>{Math.round(shareAt(mark))}%</strong><span>erram até {mark} kcal/mol</span></li>
        ))}
      </ul>
    </div>
  )
}

export interface SizeDatum { nheavy: number; n: number; mean: number; mean_abs: number }

/** Mean signed and absolute error per molecule size; bars = typical size of the error. */
export function SizeErrorChart({ data, height = 280 }: { data: SizeDatum[]; height?: number }) {
  const Content = tooltipFor<SizeDatum>((datum) => (
    <>
      <p className="chart-tooltip__title">{datum.nheavy} átomos pesados</p>
      <TooltipRow label="Moléculas" value={datum.n.toLocaleString('pt-BR')} />
      <TooltipRow label="Erro absoluto médio" value={`${formatNumber(datum.mean_abs)} kcal/mol`} color={series[0]} />
      <TooltipRow label="Erro médio (com sinal)" value={`${formatSigned(datum.mean)} kcal/mol`} color={series[1]} />
    </>
  ))
  return (
    <div className="chart" role="img" aria-label={`Erro por tamanho. ${data.map((d) => `${d.nheavy} átomos: ${formatNumber(d.mean_abs)}`).join('; ')}`}>
      <ResponsiveContainer width="100%" height={height}>
        <ComposedChart data={data} margin={margin}>
          <CartesianGrid vertical={false} stroke={chrome.grid} />
          <XAxis dataKey="nheavy" tick={axisTick} tickLine={false} axisLine={{ stroke: chrome.axis }}
            label={{ value: 'Átomos pesados (tamanho da molécula)', position: 'insideBottom', offset: -16, fill: chrome.muted, fontSize: 12 }} />
          <YAxis tick={axisTick} tickLine={false} axisLine={false} width={44} />
          <Tooltip content={Content} cursor={{ fill: 'var(--hover-wash)' }} />
          <ReferenceLine y={0} stroke={chrome.axis} />
          <Bar dataKey="mean_abs" fill={series[0]} radius={[4, 4, 0, 0]} maxBarSize={28} />
          <Line dataKey="mean" stroke={series[1]} strokeWidth={2} dot={{ r: 3, fill: series[1], stroke: 'var(--surface-1)', strokeWidth: 2 }} isAnimationActive={false} />
        </ComposedChart>
      </ResponsiveContainer>
    </div>
  )
}

export interface GroupDatum { key: string; n: number; mae: number; bias: number }

/** Horizontal bars of mean absolute error per group (e.g. source). */
export function GroupErrorChart({ data, height }: { data: GroupDatum[]; height?: number }) {
  const Content = tooltipFor<GroupDatum>((datum) => (
    <>
      <p className="chart-tooltip__title">{datum.key}</p>
      <TooltipRow label="Pares" value={datum.n.toLocaleString('pt-BR')} />
      <TooltipRow label="Erro absoluto médio" value={`${formatNumber(datum.mae)} kcal/mol`} color={series[2]} />
      <TooltipRow label="Erro médio (com sinal)" value={`${formatSigned(datum.bias)} kcal/mol`} />
    </>
  ))
  return (
    <div className="chart" role="img" aria-label={`Erro por fonte. ${data.map((d) => `${d.key}: ${formatNumber(d.mae)}`).join('; ')}`}>
      <ResponsiveContainer width="100%" height={height ?? Math.max(160, data.length * 34 + 40)}>
        <BarChart data={data} layout="vertical" margin={{ top: 4, right: 40, bottom: 20, left: 8 }}>
          <CartesianGrid horizontal={false} stroke={chrome.grid} />
          <XAxis type="number" tick={axisTick} tickLine={false} axisLine={{ stroke: chrome.axis }}
            label={{ value: 'Erro absoluto médio (kcal/mol)', position: 'insideBottom', offset: -12, fill: chrome.muted, fontSize: 12 }} />
          <YAxis type="category" dataKey="key" tick={axisTick} tickLine={false} axisLine={false} width={120} />
          <Tooltip content={Content} cursor={{ fill: 'var(--hover-wash)' }} />
          <Bar dataKey="mae" fill={series[2]} radius={[0, 4, 4, 0]} maxBarSize={18}>
            <LabelList dataKey="mae" position="right" fill="var(--text-secondary)" fontSize={11} formatter={(value: unknown) => formatNumber(Number(value), 1)} />
          </Bar>
        </BarChart>
      </ResponsiveContainer>
    </div>
  )
}


import type { ReactNode } from 'react'
import type { DatasetProfile, Histogram } from './api'
import { ColumnChart, DonutChart, HistogramChart } from './Charts'
import type { RefMark } from './Charts'
import { ChartCard } from './StoryVisuals'
import { compositionSegments, median } from './helpers'
import { formatNumber, series } from './theme'

export function SizeCard({ profile, color = series[0] }: { profile: DatasetProfile; color?: string }) {
  const typical = median(profile.nheavy)
  return (
    <ChartCard
      title="Qual o tamanho das moléculas?"
      subtitle="Quantas moléculas existem para cada número de átomos pesados (C, N, O)."
      read={<>Cada coluna é um tamanho. A molécula típica tem <strong>{typical ?? '—'} átomos pesados</strong>. Colunas altas mostram onde o conjunto tem mais exemplos — e onde um modelo terá mais dados para aprender.</>}
    >
      {Object.keys(profile.nheavy).length ? <ColumnChart counts={profile.nheavy} xLabel="Átomos pesados" color={color} /> : <p className="empty-state">Coluna nheavy ausente.</p>}
    </ChartCard>
  )
}

export function CompositionCard({ profile }: { profile: DatasetProfile }) {
  const segments = compositionSegments(profile)
  const total = segments.reduce((sum, s) => sum + s.count, 0)
  return (
    <ChartCard
      title="Quais elementos aparecem?"
      subtitle="Famílias químicas pelos elementos presentes além de C e H."
      read="Cada fatia agrupa moléculas com os mesmos tipos de átomo. Um conjunto equilibrado cobre melhor a química que o modelo precisa entender."
    >
      <DonutChart segments={segments} centerValue={total.toLocaleString('pt-BR')} centerCaption="moléculas" label="Composição por elementos" />
    </ChartCard>
  )
}

export function EnthalpyCard({ histogram, title, subtitle, color = series[0], refs, read }: {
  histogram: Histogram | undefined
  title: string
  subtitle: string
  color?: string
  refs?: RefMark[]
  read?: ReactNode
}) {
  if (!histogram || histogram.n === 0) {
    return <ChartCard title={title} subtitle={subtitle}><p className="empty-state">Sem valores numéricos nesta coluna.</p></ChartCard>
  }
  const outside = histogram.below + histogram.above
  return (
    <ChartCard
      title={title}
      subtitle={subtitle}
      read={read ?? <>Cada barra conta moléculas numa faixa de valores. Metade das moléculas está abaixo de <strong>{formatNumber(histogram.median, 1)} {histogram.unit}</strong> (mediana); 90% ficam entre {formatNumber(histogram.p05, 0)} e {formatNumber(histogram.p95, 0)}.</>}
      footer={outside > 0 ? <p className="chart-footnote">{outside.toLocaleString('pt-BR')} valores extremos (de {formatNumber(histogram.minimum, 0)} a {formatNumber(histogram.maximum, 0)}) ficam fora do eixo para não achatar o gráfico.</p> : undefined}
    >
      <HistogramChart
        histogram={histogram} color={color} xLabel={`${histogram.label} (${histogram.unit})`}
        refs={refs ?? (histogram.median !== null ? [{ value: histogram.median, label: 'mediana' }] : [])}
      />
    </ChartCard>
  )
}

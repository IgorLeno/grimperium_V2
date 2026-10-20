import type { DatasetProfile, Histogram, TrainingProgress } from './api'
import type { Segment } from './StoryVisuals'
import { series, statusColor } from './theme'

const COMPOSITION: { key: keyof DatasetProfile['composition']; label: string; description: string }[] = [
  { key: 'CH', label: 'Só C e H', description: 'Hidrocarbonetos, como o metano' },
  { key: 'CHO', label: 'C, H e O', description: 'Álcoois, ácidos, ésteres…' },
  { key: 'CHN', label: 'C, H e N', description: 'Aminas, nitrilas…' },
  { key: 'CHNO', label: 'C, H, N e O', description: 'Amidas, aminoácidos…' },
]

export function compositionSegments(profile: DatasetProfile): Segment[] {
  return COMPOSITION.map((item, index) => ({
    key: item.key, label: item.label, description: item.description,
    count: profile.composition[item.key] ?? 0, color: series[index],
  }))
}

export function findHistogram(profile: DatasetProfile, column: string): Histogram | undefined {
  return profile.histograms.find((item) => item.column === column)
}

export function median(counts: Record<string, number>): number | null {
  const entries = Object.entries(counts).map(([k, v]) => [Number(k), v] as const).sort((a, b) => a[0] - b[0])
  const total = entries.reduce((sum, [, v]) => sum + v, 0)
  let seen = 0
  for (const [value, count] of entries) {
    seen += count
    if (seen >= total / 2) return value
  }
  return null
}

/** Estados da fila de treinamento (semi_imperium.training.molecule_queue) em ordem de pipeline. */
export function queueSegments(queue: TrainingProgress): Segment[] {
  const segments: Segment[] = [
    { key: 'ok', label: 'Concluída (ok)', count: queue.ok, color: statusColor.good, icon: '✓', description: 'Protocolo 2 terminou e o mínimo foi confirmado' },
    { key: 'running', label: 'Rodando (running)', count: queue.running, color: statusColor.info, icon: '▶', description: 'Em cálculo agora' },
    { key: 'pending', label: 'Na fila (pending)', count: queue.pending, color: statusColor.neutral, icon: '…', description: 'Ainda não processada' },
    { key: 'rerun', label: 'Repetir (rerun)', count: queue.rerun, color: statusColor.warning, icon: '↻', description: 'Falhou e será tentada de novo' },
    { key: 'skip', label: 'Pulada (skip)', count: queue.skipped, color: statusColor.serious, icon: '⤼', description: 'Excedeu o máximo de tentativas' },
  ]
  return segments.filter((segment) => segment.count > 0)
}

/** Plain-language header and help text for every column the API can return. */
export const COLUMN_INFO: Record<string, { label: string; help?: string; signed?: boolean }> = {
  exp_id: { label: 'ID', help: 'Identificador na tabela experimental' },
  mol_id: { label: 'ID', help: 'Identificador da molécula' },
  name: { label: 'Nome' },
  formula: { label: 'Fórmula' },
  smiles: { label: 'SMILES', help: 'Estrutura escrita como texto' },
  nheavy: { label: 'Pesados', help: 'Átomos que não são hidrogênio' },
  charge: { label: 'Carga' },
  multiplicity: { label: 'Mult.', help: 'Multiplicidade de spin (1 = camada fechada)' },
  H298_exp: { label: 'ΔHf exp.', help: 'kcal/mol, medido' },
  uncertainty: { label: '± Incerteza', help: 'kcal/mol, quando a fonte informa' },
  source: { label: 'Fonte' },
  n_sources: { label: 'Nº fontes', help: 'Quantas fontes mediram esta molécula' },
  validation_only: { label: 'Só validação' },
  conflict_flag: { label: 'Conflito' },
  H298_cbs: { label: 'ΔHf CBS', help: 'kcal/mol, cálculo de referência' },
  H298_b3: { label: 'ΔHf B3LYP', help: 'kcal/mol, etapa intermediária do CBS' },
  cbs_b3: { label: 'CBS − B3LYP', help: 'Correção que o CBS aplica', signed: true },
  smiles_count: { label: 'Cópias', help: 'Linhas com o mesmo SMILES (conformações distintas)' },
}

/** Round tick positions (1, 2, 2.5 or 5 × 10^k) spanning [lo, hi]. */
export function niceTicks(lo: number, hi: number, target = 6): number[] {
  const span = hi - lo
  if (!(span > 0)) return [lo]
  const power = 10 ** Math.floor(Math.log10(span / target))
  // Smallest round step that keeps the tick count readable (≤ target + 2).
  const step = ([1, 2, 2.5, 5, 10].find((factor) => span / (factor * power) <= target + 2) ?? 10) * power
  const ticks: number[] = []
  for (let value = Math.ceil(lo / step) * step; value <= hi + step * 1e-6; value += step) ticks.push(Number(value.toFixed(6)) + 0)
  return ticks
}

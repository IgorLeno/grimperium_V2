import type { DatasetProfile, Histogram } from './api'
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

/** Batch states (crest_pm7.batch.enums.MoleculeStatus) in pipeline order. */
export const STATUS_INFO: { key: string; label: string; color: string; icon: string; text: string }[] = [
  { key: 'OK', label: 'Concluída', color: statusColor.good, icon: '✓', text: 'CREST e MOPAC terminaram com sucesso' },
  { key: 'Running', label: 'Rodando', color: statusColor.info, icon: '▶', text: 'Em cálculo agora' },
  { key: 'Selected', label: 'Selecionada', color: series[2], icon: '◆', text: 'Escolhida para o próximo lote' },
  { key: 'Pending', label: 'Na fila', color: statusColor.neutral, icon: '…', text: 'Ainda não processada' },
  { key: 'Rerun', label: 'Repetir', color: statusColor.warning, icon: '↻', text: 'Falhou e será tentada de novo' },
  { key: 'Skip', label: 'Pulada', color: statusColor.serious, icon: '⤼', text: 'Excedeu o máximo de tentativas' },
  { key: 'Failed', label: 'Falhou', color: statusColor.critical, icon: '✕', text: 'Erro sem nova tentativa' },
]

export function statusSegments(counts: Record<string, number>): Segment[] {
  const known = STATUS_INFO.map((info) => ({
    key: info.key, label: `${info.label} (${info.key})`, count: counts[info.key] ?? 0, color: info.color, icon: info.icon, description: info.text,
  }))
  const others = Object.entries(counts)
    .filter(([key]) => !STATUS_INFO.some((info) => info.key === key))
    .map(([key, count]) => ({ key, label: key, count, color: 'var(--status-neutral)' }))
  return [...known.filter((segment) => segment.count > 0), ...others]
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
  status: { label: 'Estado' },
  H298_pm7: { label: 'ΔHf PM7', help: 'kcal/mol, cálculo rápido' },
  delta_cbs_pm7: { label: 'Δ (CBS − PM7)', help: 'O erro que o modelo aprende', signed: true },
  cbs_quality_flag: { label: 'Qualidade CBS' },
  crest_status: { label: 'CREST' },
  mopac_status: { label: 'MOPAC' },
  crest_conformers_generated: { label: 'Conformações' },
  crest_time_s: { label: 'CREST (s)' },
  mopac_time_s: { label: 'MOPAC (s)' },
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

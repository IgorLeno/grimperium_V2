import type { TrainingStatus } from './api'

/** O coordenador grava o heartbeat a cada ≤30 s; três ciclos sem sinal é alerta. */
export const HEARTBEAT_ALERT_S = 90

export interface RunStats {
  /** Moléculas que esta rodada pretende processar (limite ou fila restante). */
  target: number
  /** Concluídas + com falha nesta rodada. */
  processed: number
  elapsedS: number
  /** Moléculas por hora (relógio de parede), `null` antes da primeira. */
  perHour: number | null
  /** Segundos até terminar, `null` enquanto não há ritmo medido. */
  etaS: number | null
}

const seconds = (from: string, to: number) => Math.max(0, (to - Date.parse(from)) / 1000)

/**
 * Ritmo e ETA da rodada atual pelo relógio de parede (inclui falhas e o
 * paralelismo real), não pelo `elapsed_s` de cada molécula.
 */
export function runStats(data: TrainingStatus, now: number): RunStats | null {
  const { job, progress } = data
  if (!job.started_at) return null
  const end = job.finished_at && !['running', 'stopping'].includes(job.state) ? Date.parse(job.finished_at) : now
  const processed = job.done + job.failed
  // A rodada vai até o limite ou até esvaziar a fila, o que vier primeiro.
  const untilEmpty = processed + progress.pending + progress.rerun + progress.running
  const target = Math.max(processed, Math.min(job.limit ?? untilEmpty, untilEmpty))
  const elapsedS = seconds(job.started_at, end)
  const perHour = processed > 0 && elapsedS > 0 ? (processed * 3600) / elapsedS : null
  const etaS = perHour ? ((target - processed) / perHour) * 3600 : null
  return { target, processed, elapsedS, perHour, etaS }
}

export function heartbeatAge(data: TrainingStatus, now: number): number | null {
  return data.job.heartbeat_at ? seconds(data.job.heartbeat_at, now) : null
}

export function elapsedSince(iso: string | null, now: number): number | null {
  return iso ? seconds(iso, now) : null
}

/** `42 s`, `3 min 05 s`, `2 h 07 min`. */
export function formatDuration(totalS: number): string {
  const s = Math.round(totalS)
  if (s < 60) return `${s} s`
  const m = Math.floor(s / 60)
  if (m < 60) return `${m} min ${String(s % 60).padStart(2, '0')} s`
  const h = Math.floor(m / 60)
  if (h < 48) return `${h} h ${String(m % 60).padStart(2, '0')} min`
  return `${Math.floor(h / 24)} d ${h % 24} h`
}

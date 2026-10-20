import { expect, test } from 'vitest'
import type { TrainingStatus } from './api'
import { formatDuration, heartbeatAge, runStats } from './trainingRun'

const T0 = Date.parse('2026-10-20T10:00:00Z')

function status(job: Partial<TrainingStatus['job']>, queued = 100): TrainingStatus {
  return {
    progress: { total: 1000, ok: 0, pending: queued, rerun: 0, running: 0, skipped: 0 },
    job: { state: 'running', pid: 1, started_at: '2026-10-20T10:00:00Z', heartbeat_at: '2026-10-20T10:00:00Z', finished_at: null, workers: 4, limit: 50, done: 0, failed: 0, message: '', stop_requested: false, ...job },
    geometries: 0,
    failures: [],
    activity: { in_flight: [], recent: [], mean_elapsed_s: null },
    cpu_count: 16,
  }
}

test('sem rodada não há estatística', () => {
  expect(runStats(status({ started_at: null }), T0)).toBeNull()
})

test('antes da primeira molécula não inventa ritmo nem ETA', () => {
  const stats = runStats(status({}), T0 + 600_000)!
  expect(stats).toMatchObject({ target: 50, processed: 0, perHour: null, etaS: null })
  expect(stats.elapsedS).toBe(600)
})

test('ritmo e ETA pelo relógio de parede, falhas contam como processadas', () => {
  const stats = runStats(status({ done: 8, failed: 2 }), T0 + 3_600_000)!
  expect(stats.perHour).toBe(10)
  expect(stats.etaS).toBe(4 * 3600) // 40 restantes a 10/h
})

test('sem limite o alvo é o processado mais a fila restante', () => {
  const stats = runStats(status({ limit: null, done: 5 }, 95), T0 + 60_000)!
  expect(stats.target).toBe(100)
})

test('rodada encerrada para o relógio no finished_at', () => {
  const stats = runStats(status({ state: 'finished', done: 50, finished_at: '2026-10-20T12:00:00Z' }), T0 + 86_400_000)!
  expect(stats.elapsedS).toBe(7200)
  expect(stats.etaS).toBe(0)
})

test('idade do heartbeat', () => {
  expect(heartbeatAge(status({}), T0 + 12_000)).toBe(12)
  expect(heartbeatAge(status({ heartbeat_at: null }), T0)).toBeNull()
})

test('formata durações', () => {
  expect(formatDuration(42)).toBe('42 s')
  expect(formatDuration(185)).toBe('3 min 05 s')
  expect(formatDuration(7620)).toBe('2 h 07 min')
  expect(formatDuration(3 * 86400 + 3600)).toBe('3 d 1 h')
})

test('limite maior que a fila restante para na fila', () => {
  expect(runStats(status({ limit: 500, done: 2 }, 8), T0 + 60_000)!.target).toBe(10)
})

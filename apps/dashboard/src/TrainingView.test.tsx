import { cleanup, fireEvent, render, screen, waitFor } from '@testing-library/react'
import { afterEach, expect, test, vi } from 'vitest'
import { TrainingView } from './TrainingView'
import type { TrainingStatus } from './api'

const progress = {
  total: 29568, protocol2: 10, protocol1: 8232, ok_unmarked: 0, pending: 21000,
  selected: 0, rerun: 3, running: 1, skipped: 29, migrated: true,
}
const idle: TrainingStatus['job'] = { state: 'idle', pid: null, started_at: null, heartbeat_at: null, finished_at: null, workers: 0, limit: null, done: 0, failed: 0, message: '', stop_requested: false }

function status(overrides: Partial<TrainingStatus> = {}): TrainingStatus {
  return { progress, job: idle, geometries: 10, failures: [], ...overrides }
}

function json(data: unknown, code = 200): Response {
  return new Response(JSON.stringify(data), { status: code, headers: { 'Content-Type': 'application/json' } })
}

function mock(current: () => TrainingStatus) {
  const calls: { path: string; method: string; body: string | null }[] = []
  vi.stubGlobal('fetch', vi.fn((input: RequestInfo | URL, init?: RequestInit) => {
    const path = String(input)
    calls.push({ path, method: init?.method ?? 'GET', body: typeof init?.body === 'string' ? init.body : null })
    return Promise.resolve(json(current()))
  }))
  return calls
}

afterEach(() => { cleanup(); vi.unstubAllGlobals() })

test('mostra as contagens com denominador explícito', async () => {
  mock(() => status())
  render(<TrainingView />)
  const bar = await screen.findByRole('img', { name: /Andamento do treinamento/ })
  expect(bar).toHaveAccessibleName(/Protocolo 2: 10 de 29\.568/)
  expect(screen.getAllByText(/8\.232 preliminares/).length).toBeGreaterThan(0)
  expect(screen.getAllByText(/21\.003 na fila/).length).toBeGreaterThan(0)
})

test('inicia com workers e limite e mostra o estado do job', async () => {
  let current = status()
  const calls = mock(() => current)
  render(<TrainingView />)
  await screen.findByRole('img', { name: /Andamento do treinamento/ })
  fireEvent.change(screen.getByLabelText('Moléculas em paralelo'), { target: { value: '2' } })
  fireEvent.change(screen.getByLabelText('Limite de moléculas (opcional)'), { target: { value: '5' } })
  current = status({ job: { ...idle, state: 'running', workers: 2, limit: 5, pid: 123 } })
  fireEvent.click(screen.getByRole('button', { name: 'Iniciar treinamento' }))
  await waitFor(() => expect(calls.some((c) => c.path === '/api/training/start')).toBe(true))
  const start = calls.find((c) => c.path === '/api/training/start')!
  expect(start.method).toBe('POST')
  expect(JSON.parse(start.body ?? '{}')).toEqual({ workers: 2, limit: 5 })
  expect(await screen.findByText(/Em execução com 2 em paralelo/)).toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Parar' })).toBeEnabled()
  expect(screen.getByRole('button', { name: 'Iniciar treinamento' })).toBeDisabled()
})

test('bloqueia o início quando o CSV não foi migrado e mostra o comando', async () => {
  mock(() => status({ progress: { ...progress, migrated: false, ok_unmarked: 8232 } }))
  render(<TrainingView />)
  await screen.findByRole('img', { name: /Andamento do treinamento/ })
  expect(screen.getByRole('button', { name: 'Iniciar treinamento' })).toBeDisabled()
  expect(screen.getByText(/semi_imperium\.training migrate --apply/)).toBeInTheDocument()
})

test('lista as falhas recentes com o motivo', async () => {
  mock(() => status({ failures: [{ mol_id: 'mol_00042', reason: 'saddle: imaginary mode', at: '2026-10-12T10:00:00+00:00' }] }))
  render(<TrainingView />)
  expect(await screen.findByText('mol_00042')).toBeInTheDocument()
  expect(screen.getByText(/saddle: imaginary mode/)).toBeInTheDocument()
})

test('mostra o erro da API ao iniciar', async () => {
  vi.stubGlobal('fetch', vi.fn((input: RequestInfo | URL) =>
    Promise.resolve(String(input) === '/api/training/start' ? json({ detail: 'Já existe um treinamento em execução.' }, 409) : json(status()))))
  render(<TrainingView />)
  await screen.findByRole('img', { name: /Andamento do treinamento/ })
  fireEvent.click(screen.getByRole('button', { name: 'Iniciar treinamento' }))
  expect(await screen.findByRole('alert')).toHaveTextContent('Já existe um treinamento em execução.')
})

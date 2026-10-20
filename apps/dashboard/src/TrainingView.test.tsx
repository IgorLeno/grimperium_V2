import { cleanup, fireEvent, render, screen, waitFor } from '@testing-library/react'
import { afterEach, expect, test, vi } from 'vitest'
import { TrainingView } from './TrainingView'
import type { TrainingStatus } from './api'

const progress = { total: 27335, ok: 10, pending: 21000, rerun: 3, running: 1, skipped: 29 }
const idle: TrainingStatus['job'] = { state: 'idle', pid: null, started_at: null, heartbeat_at: null, finished_at: null, workers: 0, limit: null, done: 0, failed: 0, message: '', stop_requested: false }

function status(overrides: Partial<TrainingStatus> = {}): TrainingStatus {
  return { progress, job: idle, geometries: 10, failures: [], activity: { in_flight: [], recent: [], mean_elapsed_s: null }, cpu_count: 16, ...overrides }
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
  expect(bar).toHaveAccessibleName(/Concluídas: 10 de 27\.335/)
  expect(screen.getAllByText(/21\.003 na fila/).length).toBeGreaterThan(0)
  expect(screen.getByText(/10 de 27\.335 moléculas concluídas/)).toBeInTheDocument()
  expect(screen.queryByText(/preliminar|protocolo 1|migra/i)).not.toBeInTheDocument()
})

test('inicia com workers, threads e limite e mostra o estado do job', async () => {
  let current = status()
  const calls = mock(() => current)
  render(<TrainingView />)
  await screen.findByRole('img', { name: /Andamento do treinamento/ })
  fireEvent.change(screen.getByLabelText('Moléculas em paralelo'), { target: { value: '2' } })
  fireEvent.change(screen.getByLabelText('Threads por molécula'), { target: { value: '6' } })
  expect(screen.getByText(/2 × 6 =/)).toHaveTextContent('CREST: 2 × 6 = 12 de 16 núcleos')
  fireEvent.change(screen.getByLabelText('Limite de moléculas (opcional)'), { target: { value: '5' } })
  current = status({ job: { ...idle, state: 'running', workers: 2, limit: 5, pid: 123 } })
  fireEvent.click(screen.getByRole('button', { name: 'Iniciar cálculo' }))
  await waitFor(() => expect(calls.some((c) => c.path === '/api/training/start')).toBe(true))
  const start = calls.find((c) => c.path === '/api/training/start')!
  expect(start.method).toBe('POST')
  expect(JSON.parse(start.body ?? '{}')).toEqual({ workers: 2, crest_threads: 6, limit: 5 })
  expect(await screen.findByText(/Em execução com 2 em paralelo/)).toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Parar' })).toBeEnabled()
  expect(screen.getByRole('button', { name: 'Iniciar cálculo' })).toBeDisabled()
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
  fireEvent.click(screen.getByRole('button', { name: 'Iniciar cálculo' }))
  expect(await screen.findByRole('alert')).toHaveTextContent('Já existe um treinamento em execução.')
})

test('avisa quando workers × threads passa dos núcleos', async () => {
  mock(() => status())
  render(<TrainingView />)
  await screen.findByRole('img', { name: /Andamento do treinamento/ })
  fireEvent.change(screen.getByLabelText('Moléculas em paralelo'), { target: { value: '4' } })
  fireEvent.change(screen.getByLabelText('Threads por molécula'), { target: { value: '8' } })
  expect(screen.getByText(/acima dos núcleos disponíveis/)).toHaveTextContent('32 de 16 núcleos')
})

test('rodada ativa: barra da rodada, ritmo, cartões por worker e feed', async () => {
  const started = new Date(Date.now() - 3_600_000).toISOString()
  const job = { ...idle, state: 'running' as const, pid: 1, workers: 2, limit: 50, done: 8, failed: 2, started_at: started, heartbeat_at: new Date().toISOString() }
  mock(() => status({
    job,
    activity: {
      in_flight: [{ mol_id: 'mol_00007', smiles: 'CCO', nheavy: 3, started_at: new Date(Date.now() - 30_000).toISOString() }],
      recent: [{ mol_id: 'mol_00003', smiles: 'CC', h298_pm7: -20.04, delta_kcal_mol: 0.5, elapsed_s: 125, conformers_evaluated: 4, finished_at: new Date().toISOString() }],
      mean_elapsed_s: 300,
    },
    failures: [{ mol_id: 'mol_00042', reason: 'saddle: imaginary mode', at: started }],
  }))
  render(<TrainingView />)
  const bar = await screen.findByRole('progressbar', { name: /Rodada: 10 de 50/ })
  expect(bar).toHaveAttribute('aria-valuenow', '10')
  expect(screen.getByText(/10 mol\/h/)).toBeInTheDocument()
  expect(screen.getByText(/~4 h 00 min/)).toBeInTheDocument()
  expect(screen.getByText(/Coordenador ativo/)).toBeInTheDocument()
  expect(screen.getByText('mol_00007')).toBeInTheDocument()
  expect(screen.getByText(/Aguardando a próxima molécula/)).toBeInTheDocument()
  expect(screen.getByText(/ΔHf PM7 -20,0 kcal\/mol/)).toBeInTheDocument()
  expect(screen.getByText('mol_00042')).toBeInTheDocument()
})

test('heartbeat parado vira alerta', async () => {
  const old = new Date(Date.now() - 5 * 60_000).toISOString()
  mock(() => status({ job: { ...idle, state: 'running', pid: 1, workers: 1, started_at: old, heartbeat_at: old } }))
  render(<TrainingView />)
  expect(await screen.findByText(/Sem sinal do coordenador há 5 min/)).toBeInTheDocument()
})

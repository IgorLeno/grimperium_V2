import { useCallback, useEffect, useState } from 'react'
import { getJson } from './api'

type ApiState<T> =
  | { status: 'loading' }
  | { status: 'error'; message: string }
  | { status: 'ready'; data: T }

interface Options {
  /** Atualiza a cada `pollMs` ms sem voltar ao estado "carregando". */
  pollMs?: number
}

export function useApi<T>(path: string, options: Options = {}): { state: ApiState<T>; reload: () => void; refresh: () => void } {
  const { pollMs } = options
  const [state, setState] = useState<ApiState<T>>({ status: 'loading' })
  const [revision, setRevision] = useState(0)

  useEffect(() => {
    const controller = new AbortController()
    const load = () => getJson<T>(path, controller.signal).then(
      (data) => setState({ status: 'ready', data }),
      (error: unknown) => {
        if (!controller.signal.aborted) {
          setState({ status: 'error', message: String(error instanceof Error ? error.message : error) })
        }
      },
    )
    void load()
    const timer = pollMs ? window.setInterval(() => { void load() }, pollMs) : undefined
    return () => {
      controller.abort()
      if (timer !== undefined) window.clearInterval(timer)
    }
  }, [path, revision, pollMs])

  const reload = useCallback(() => {
    setState({ status: 'loading' })
    setRevision((current) => current + 1)
  }, [])

  // Busca de novo sem passar por "carregando": a tela não pisca.
  const refresh = useCallback(() => setRevision((current) => current + 1), [])

  return { state, reload, refresh }
}

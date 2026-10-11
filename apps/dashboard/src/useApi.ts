import { useCallback, useEffect, useState } from 'react'
import { getJson } from './api'

type ApiState<T> =
  | { status: 'loading' }
  | { status: 'error'; message: string }
  | { status: 'ready'; data: T }

export function useApi<T>(path: string): { state: ApiState<T>; reload: () => void } {
  const [state, setState] = useState<ApiState<T>>({ status: 'loading' })
  const [revision, setRevision] = useState(0)

  useEffect(() => {
    const controller = new AbortController()
    getJson<T>(path, controller.signal).then(
      (data) => setState({ status: 'ready', data }),
      (error: unknown) => {
        if (!controller.signal.aborted) {
          setState({ status: 'error', message: String(error instanceof Error ? error.message : error) })
        }
      },
    )
    return () => controller.abort()
  }, [path, revision])

  const reload = useCallback(() => {
    setState({ status: 'loading' })
    setRevision((current) => current + 1)
  }, [])

  return { state, reload }
}

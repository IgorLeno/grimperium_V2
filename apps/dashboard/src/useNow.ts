import { useEffect, useState } from 'react'

/** Relógio de 1 s para cronômetros e idade do heartbeat entre um poll e outro. */
export function useNow(active: boolean): number {
  const [now, setNow] = useState(() => Date.now())
  useEffect(() => {
    if (!active) return undefined
    const id = window.setInterval(() => setNow(Date.now()), 1000)
    return () => window.clearInterval(id)
  }, [active])
  return now
}

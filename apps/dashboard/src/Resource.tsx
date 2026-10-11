import type { ReactNode } from 'react'

export function Section({ title, children }: { title: string; children: ReactNode }) {
  return (
    <section className="panel" aria-label={title}>
      <h2>{title}</h2>
      {children}
    </section>
  )
}

export function ErrorNotice({ message, retry }: { message: string; retry: () => void }) {
  return (
    <div role="alert" className="notice error">
      <p>{message}</p>
      <button type="button" onClick={retry}>Tentar novamente</button>
    </div>
  )
}

export function CountTable({
  counts,
  label,
}: {
  counts: Record<string, number>
  label: string
}) {
  const entries = Object.entries(counts).sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
  if (entries.length === 0) return <p>Nenhuma contagem disponível.</p>
  return (
    <table>
      <thead><tr><th scope="col">{label}</th><th scope="col">Quantidade</th></tr></thead>
      <tbody>
        {entries.map(([name, count]) => (
          <tr key={name}><th scope="row">{name}</th><td>{count.toLocaleString('pt-BR')}</td></tr>
        ))}
      </tbody>
    </table>
  )
}

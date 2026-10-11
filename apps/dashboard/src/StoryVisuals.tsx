import type { ReactNode } from 'react'

export function MoleculeIllustration() {
  return (
    <svg className="molecule-art" viewBox="0 0 420 270" fill="none" aria-hidden="true">
      <circle cx="210" cy="135" r="112" stroke="currentColor" strokeOpacity=".18" />
      <circle cx="210" cy="135" r="79" stroke="currentColor" strokeOpacity=".18" strokeDasharray="4 8" />
      <path d="M101 174 161 99 231 124 291 76 343 144 264 206 191 178 101 174Z" stroke="currentColor" strokeOpacity=".55" strokeWidth="2" />
      <path d="m161 99 30 79m40-54 33 82m27-130-27 130" stroke="currentColor" strokeOpacity=".28" strokeWidth="2" />
      <circle cx="101" cy="174" r="10" fill="#85d9cb" />
      <circle cx="161" cy="99" r="14" fill="#e9bd7b" />
      <circle cx="231" cy="124" r="19" fill="#f7f2e8" />
      <circle cx="291" cy="76" r="11" fill="#85d9cb" />
      <circle cx="343" cy="144" r="15" fill="#e9bd7b" />
      <circle cx="264" cy="206" r="12" fill="#f7f2e8" />
      <circle cx="191" cy="178" r="8" fill="#85d9cb" />
      <circle cx="210" cy="135" r="3" fill="#f7f2e8" />
    </svg>
  )
}

export function MetricCard({ label, value, note, accent = false }: {
  label: string
  value: string
  note: string
  accent?: boolean
}) {
  return (
    <div className={`metric-card${accent ? ' metric-card--accent' : ''}`}>
      <p className="metric-card__label">{label}</p>
      <strong className="metric-card__value">{value}</strong>
      <p className="metric-card__note">{note}</p>
    </div>
  )
}

export function RankedBars({ counts, total, label, limit = 6 }: {
  counts: Record<string, number>
  total: number
  label: string
  limit?: number
}) {
  const entries = Object.entries(counts).sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
  if (entries.length === 0) return <p className="empty-state">Nenhuma contagem disponível.</p>
  const top = entries.slice(0, limit)
  const remaining = entries.slice(limit).reduce((sum, [, count]) => sum + count, 0)
  const bars = remaining > 0 ? [...top, ['Outras fontes', remaining] as [string, number]] : top
  const max = Math.max(...bars.map(([, count]) => count), 1)
  return (
    <div className="ranked-bars" role="img" aria-label={`${label}. ${entries.map(([name, count]) => `${name}: ${count}`).join('; ')}`}>
      {bars.map(([name, count], index) => (
        <div className="ranked-bars__row" key={name}>
          <span className="ranked-bars__name" title={name}>{name}</span>
          <span className="ranked-bars__track"><span className={index === 0 ? 'ranked-bars__fill ranked-bars__fill--lead' : 'ranked-bars__fill'} style={{ width: `${(count / max) * 100}%` }} /></span>
          <span className="ranked-bars__count">{count.toLocaleString('pt-BR')}</span>
          <span className="ranked-bars__share">{total > 0 ? `${Math.round((count / total) * 100)}%` : '—'}</span>
        </div>
      ))}
    </div>
  )
}

export function StoryNote({ number, title, children }: { number: string; title: string; children: ReactNode }) {
  return (
    <div className="story-note">
      <span className="story-note__number">{number}</span>
      <div><h3>{title}</h3><p>{children}</p></div>
    </div>
  )
}

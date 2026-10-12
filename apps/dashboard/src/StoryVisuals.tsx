import type { ReactNode } from 'react'

export function MoleculeIllustration() {
  return (
    <svg className="molecule-art" viewBox="0 0 420 300" fill="none" aria-hidden="true">
      <defs>
        <radialGradient id="atom-a" cx="35%" cy="30%" r="70%">
          <stop offset="0" stopColor="#d9fff4" />
          <stop offset="1" stopColor="#3fc6a8" />
        </radialGradient>
        <radialGradient id="atom-b" cx="35%" cy="30%" r="70%">
          <stop offset="0" stopColor="#fff1dc" />
          <stop offset="1" stopColor="#f0a35e" />
        </radialGradient>
        <radialGradient id="atom-c" cx="35%" cy="30%" r="70%">
          <stop offset="0" stopColor="#ffffff" />
          <stop offset="1" stopColor="#b9cfe0" />
        </radialGradient>
      </defs>
      <g className="molecule-art__orbits">
        <ellipse cx="210" cy="150" rx="168" ry="62" stroke="currentColor" strokeOpacity=".16" transform="rotate(-18 210 150)" />
        <ellipse cx="210" cy="150" rx="168" ry="62" stroke="currentColor" strokeOpacity=".12" transform="rotate(32 210 150)" />
        <circle cx="210" cy="150" r="118" stroke="currentColor" strokeOpacity=".1" strokeDasharray="3 9" />
      </g>
      <g className="molecule-art__bonds" stroke="currentColor" strokeOpacity=".5" strokeWidth="3" strokeLinecap="round">
        <path d="M108 190 160 112 236 136 296 84 348 156 270 216 196 194 108 190" />
        <path d="m160 112 36 82m40-58 34 80" strokeOpacity=".28" />
      </g>
      <g className="molecule-art__atoms">
        <circle cx="108" cy="190" r="12" fill="url(#atom-a)" />
        <circle cx="160" cy="112" r="17" fill="url(#atom-b)" />
        <circle cx="236" cy="136" r="23" fill="url(#atom-c)" />
        <circle cx="296" cy="84" r="13" fill="url(#atom-a)" />
        <circle cx="348" cy="156" r="18" fill="url(#atom-b)" />
        <circle cx="270" cy="216" r="14" fill="url(#atom-c)" />
        <circle cx="196" cy="194" r="10" fill="url(#atom-a)" />
      </g>
    </svg>
  )
}

export type Tone = 'default' | 'accent' | 'good' | 'warning' | 'critical'

export function MetricCard({ label, value, note, icon, unit, accent = false, tone, hint }: {
  label: string
  value: string
  note: string
  icon?: ReactNode
  unit?: string
  accent?: boolean
  tone?: Tone
  /** One sentence for a lay reader: what the number means. */
  hint?: string
}) {
  const resolved = tone ?? (accent ? 'accent' : 'default')
  return (
    <div className={`metric-card metric-card--${resolved}${accent ? ' metric-card--accent' : ''}`}>
      <div className="metric-card__top">
        <p className="metric-card__label">{label}</p>
        {icon && <span className="metric-card__icon" aria-hidden="true">{icon}</span>}
      </div>
      <strong className="metric-card__value">{value}{unit && <small>{unit}</small>}</strong>
      <p className="metric-card__note">{note}</p>
      {hint && <p className="metric-card__hint"><span aria-hidden="true">💡</span> {hint}</p>}
    </div>
  )
}

export function ChartCard({ title, subtitle, read, children, wide = false, footer }: {
  title: string
  subtitle?: string
  /** "Como ler" — a plain-language guide to the chart. */
  read?: ReactNode
  children: ReactNode
  wide?: boolean
  footer?: ReactNode
}) {
  return (
    <section className={`chart-card${wide ? ' chart-card--wide' : ''}`} aria-label={title}>
      <header className="chart-card__header">
        <h3>{title}</h3>
        {subtitle && <p>{subtitle}</p>}
      </header>
      <div className="chart-card__body">{children}</div>
      {read && <div className="how-to-read"><span>Como ler</span><p>{read}</p></div>}
      {footer}
    </section>
  )
}

export function RingGauge({ share, label, caption, color = 'var(--series-1)' }: {
  /** 0–1. */
  share: number | null
  label: string
  caption: string
  color?: string
}) {
  const radius = 52
  const circumference = 2 * Math.PI * radius
  const value = share === null ? 0 : Math.max(0, Math.min(1, share))
  return (
    <div className="ring-gauge" role="img" aria-label={`${label}: ${share === null ? 'sem dados' : `${Math.round(value * 100)}%`}`}>
      <svg viewBox="0 0 128 128" aria-hidden="true">
        <circle cx="64" cy="64" r={radius} className="ring-gauge__track" />
        <circle
          cx="64" cy="64" r={radius} className="ring-gauge__value" stroke={color}
          strokeDasharray={`${value * circumference} ${circumference}`} transform="rotate(-90 64 64)"
        />
      </svg>
      <div className="ring-gauge__center">
        <strong>{share === null ? '—' : `${Math.round(value * 100)}%`}</strong>
        <span>{caption}</span>
      </div>
    </div>
  )
}

export interface Segment {
  key: string
  label: string
  count: number
  color: string
  description?: string
  icon?: string
}

export function ProportionBar({ segments, total, label, showLegend = true }: {
  segments: Segment[]
  total?: number
  label: string
  showLegend?: boolean
}) {
  const sum = total ?? segments.reduce((acc, segment) => acc + segment.count, 0)
  const visible = segments.filter((segment) => segment.count > 0)
  return (
    <div className="proportion">
      <div
        className="proportion__bar" role="img"
        aria-label={`${label}. ${segments.map((segment) => `${segment.label}: ${segment.count}`).join('; ')}`}
      >
        {visible.map((segment) => (
          <span
            key={segment.key}
            title={`${segment.label}: ${segment.count.toLocaleString('pt-BR')} (${sum ? Math.round((segment.count / sum) * 100) : 0}%)`}
            style={{ flexGrow: segment.count, backgroundColor: segment.color }}
          >
            {sum > 0 && segment.count / sum >= 0.08 && <em>{Math.round((segment.count / sum) * 100)}%</em>}
          </span>
        ))}
      </div>
      {showLegend && (
        <ul className="proportion__legend">
          {segments.map((segment) => (
            <li key={segment.key}>
              <i style={{ backgroundColor: segment.color }} aria-hidden="true" />
              <span className="proportion__name">{segment.icon && <b aria-hidden="true">{segment.icon} </b>}{segment.label}</span>
              <strong>{segment.count.toLocaleString('pt-BR')}</strong>
              <span className="proportion__share">{sum ? `${Math.round((segment.count / sum) * 100)}%` : '—'}</span>
              {segment.description && <small>{segment.description}</small>}
            </li>
          ))}
        </ul>
      )}
    </div>
  )
}

export function RankedBars({ counts, total, label, limit = 6, color }: {
  counts: Record<string, number>
  total: number
  label: string
  limit?: number
  color?: string
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
          <span className="ranked-bars__track">
            <span
              className={index === 0 ? 'ranked-bars__fill ranked-bars__fill--lead' : 'ranked-bars__fill'}
              style={{ width: `${(count / max) * 100}%`, ...(color ? { backgroundColor: color } : {}) }}
            />
          </span>
          <span className="ranked-bars__count">{count.toLocaleString('pt-BR')}</span>
          <span className="ranked-bars__share">{total > 0 ? `${Math.round((count / total) * 100)}%` : '—'}</span>
        </div>
      ))}
    </div>
  )
}

export function StoryNote({ number, title, children, icon }: { number: string; title: string; children: ReactNode; icon?: string }) {
  return (
    <div className="story-note">
      <span className="story-note__number" aria-hidden={icon ? 'true' : undefined}>{icon ?? number}</span>
      <div><h3>{title}</h3><p>{children}</p></div>
    </div>
  )
}

const pipelineSteps = [
  { icon: '⌬', title: 'Molécula', tag: 'SMILES', text: 'A estrutura é escrita como texto (SMILES), por exemplo CCO para o etanol.' },
  { icon: '⟳', title: 'Conformações', tag: 'CREST', text: 'A molécula pode se dobrar de vários jeitos; buscamos as formas de menor energia.' },
  { icon: '⚡', title: 'Cálculo rápido', tag: 'PM7 · MOPAC', text: 'Um método barato estima a energia em segundos — rápido, porém com erro.' },
  { icon: '◎', title: 'Referência precisa', tag: 'CBS', text: 'Um método caro e rigoroso serve de gabarito para medir o erro do PM7.' },
  { icon: '✦', title: 'Correção por IA', tag: 'Δ-learning', text: 'O modelo aprende o erro (Δ) e o soma ao PM7: rapidez com precisão.' },
]

export function PipelineDiagram() {
  return (
    <ol className="pipeline" aria-label="Etapas do pipeline Grimperium">
      {pipelineSteps.map((step, index) => (
        <li key={step.title} className="pipeline__step" style={{ animationDelay: `${index * 70}ms` }}>
          <span className="pipeline__icon" aria-hidden="true">{step.icon}</span>
          <span className="pipeline__index">Passo {index + 1}</span>
          <strong>{step.title}</strong>
          <span className="pipeline__tag">{step.tag}</span>
          <p>{step.text}</p>
        </li>
      ))}
    </ol>
  )
}

/** Conceptual figure of delta-learning; the numbers are illustrative, not data. */
export function DeltaExplainer() {
  const rows = [
    { label: 'PM7 (rápido)', value: 62, display: '−48,0', className: 'delta-bar--pm7' },
    { label: 'Correção Δ aprendida', value: 18, offset: 62, display: '−6,5', className: 'delta-bar--delta' },
    { label: 'Previsão final ≈ CBS', value: 80, display: '−54,5', className: 'delta-bar--final' },
  ]
  return (
    <figure className="delta-explainer">
      <div className="delta-explainer__rows" role="img" aria-label="Ilustração: o valor PM7 somado à correção aprendida aproxima o valor CBS">
        {rows.map((row) => (
          <div key={row.label} className="delta-row">
            <span className="delta-row__label">{row.label}</span>
            <span className="delta-row__track">
              <span className={`delta-bar ${row.className}`} style={{ marginLeft: `${row.offset ?? 0}%`, width: `${row.value}%` }} />
            </span>
            <span className="delta-row__value">{row.display}</span>
          </div>
        ))}
      </div>
      <div className="delta-equation" aria-hidden="true">
        <span className="chip chip--pm7">PM7</span><b>+</b><span className="chip chip--delta">Δ previsto</span><b>=</b><span className="chip chip--final">≈ CBS</span>
      </div>
      <figcaption>Ilustração conceitual com valores fictícios (kcal/mol). O modelo não prevê a energia do zero: aprende quanto o PM7 costuma errar.</figcaption>
    </figure>
  )
}

export interface GlossaryTerm {
  term: string
  definition: string
  example?: string
}

export function Glossary({ terms, title = 'Glossário rápido' }: { terms: GlossaryTerm[]; title?: string }) {
  return (
    <section className="glossary" aria-label={title}>
      <div className="glossary__heading">
        <span className="eyebrow eyebrow--dark">Para quem está chegando</span>
        <h2>{title}</h2>
      </div>
      <dl className="glossary__grid">
        {terms.map((item) => (
          <div key={item.term} className="glossary__item">
            <dt>{item.term}</dt>
            <dd>{item.definition}{item.example && <span className="glossary__example">{item.example}</span>}</dd>
          </div>
        ))}
      </dl>
    </section>
  )
}

export function PageHero({ eyebrow, title, accent, lead, onRefresh, art, variant = 'default' }: {
  eyebrow: string
  title: string
  accent: string
  lead: ReactNode
  onRefresh?: () => void
  art?: ReactNode
  variant?: 'default' | 'exp' | 'cbs' | 'pm7' | 'control'
}) {
  return (
    <div className={`hero hero--${variant}`}>
      <div className="hero__copy">
        <span className="eyebrow">{eyebrow}</span>
        <h1>{title}<br /><em>{accent}</em></h1>
        <p>{lead}</p>
        {onRefresh && <button type="button" className="refresh-button" onClick={onRefresh}>↻ <span>Atualizar dados</span></button>}
      </div>
      <div className="hero__art">{art ?? <MoleculeIllustration />}</div>
    </div>
  )
}

export function SectionHeading({ eyebrow, title, children, badge }: {
  eyebrow: string
  title: string
  children?: ReactNode
  badge?: ReactNode
}) {
  return (
    <div className="section-heading">
      <div>
        <span className="eyebrow eyebrow--dark">{eyebrow}</span>
        <h2>{title}</h2>
        {children && <p>{children}</p>}
      </div>
      {badge && <span className="badge">{badge}</span>}
    </div>
  )
}

export interface Fact {
  icon: string
  title: string
  text: ReactNode
}

/** "What is this table?" — origin, size, purpose, in plain words. */
export function Explainer({ title, intro, facts, analogy }: {
  title: string
  intro: ReactNode
  facts: Fact[]
  analogy?: ReactNode
}) {
  return (
    <section className="explainer" aria-label={title}>
      <div className="explainer__intro">
        <span className="eyebrow eyebrow--dark">Entenda a tabela</span>
        <h2>{title}</h2>
        <p>{intro}</p>
        {analogy && <p className="explainer__analogy"><span aria-hidden="true">🧭</span> {analogy}</p>}
      </div>
      <ul className="explainer__facts">
        {facts.map((fact) => (
          <li key={fact.title}>
            <span className="explainer__icon" aria-hidden="true">{fact.icon}</span>
            <div><strong>{fact.title}</strong><p>{fact.text}</p></div>
          </li>
        ))}
      </ul>
    </section>
  )
}

export function Callout({ label = 'Destaque', children, tone = 'info' }: { label?: string; children: ReactNode; tone?: 'info' | 'warning' }) {
  return <div className={`callout callout--${tone}`}><span>{label}</span><p>{children}</p></div>
}

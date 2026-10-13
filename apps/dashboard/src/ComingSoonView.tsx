import { PageHero, SectionHeading } from './StoryVisuals'

/** Placeholder honesto de uma rotina cujo projeto ainda não foi construído. */
export function ComingSoonView({ eyebrow, title, accent, lead, depends }: {
  eyebrow: string; title: string; accent: string; lead: string; depends: string
}) {
  return (
    <div className="view">
      <PageHero variant="control" eyebrow={eyebrow} title={title} accent={accent} lead={lead} />
      <section className="panel" aria-label="Rotina em construção">
        <SectionHeading eyebrow="Status" title="Em breve" badge="Em construção">{depends}</SectionHeading>
      </section>
    </div>
  )
}

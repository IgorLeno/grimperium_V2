import { PageHero, SectionHeading } from './StoryVisuals'

/**
 * Resumos do hub: compilado das propriedades obtidas pelas rotinas.
 * Ainda não há back-end de jobs nem resultados, então a tela não mostra
 * números: apenas explica o que vai aparecer aqui.
 */
export function ResultsView() {
  return (
    <div className="view">
      <PageHero
        variant="control" eyebrow="Resumos · Resultados do hub" title="O que as rotinas" accent="já calcularam."
        lead="Aqui ficará o compilado das propriedades obtidas pelas rotinas criadas no Delta Chemistry, reunidas em um só lugar."
      />
      <section className="panel" aria-label="Resultados das rotinas">
        <SectionHeading eyebrow="Resultados" title="Nenhum resultado ainda" badge={<span className="badge">Sem cálculos</span>}>
          Nenhuma rotina foi executada, então não há valores para resumir.
        </SectionHeading>
        <div className="results-empty">
          <span className="results-empty__icon" aria-hidden="true">Δ</span>
          <div>
            <strong>Quando uma rotina rodar, os valores aparecem aqui.</strong>
            <p>O resumo reunirá, por propriedade, as moléculas calculadas, o valor rápido e o valor corrigido pelo modelo. Acompanhe o cálculo das moléculas em Delta Entalpia → Calcular PM7, na barra lateral.</p>
          </div>
        </div>
      </section>
    </div>
  )
}

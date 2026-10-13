import { MODULES } from './modules'
import { Callout, DeltaExplainer, Glossary, PageHero, SectionHeading, StoryNote } from './StoryVisuals'

const GLOSSARY = [
  { term: 'Método rápido', definition: 'Cálculo barato, de segundos por molécula, mas com erro maior (por exemplo, o PM7).' },
  { term: 'Método de referência', definition: 'Cálculo mais preciso e bem mais caro (por exemplo, o CBS). Serve de gabarito para medir o erro do método rápido.' },
  { term: 'Δ-learning', definition: 'A IA aprende só a diferença entre a referência e o método rápido (Δ), e não o valor inteiro. Precisa de menos dados.' },
  { term: 'Módulo', definition: 'Um método rápido corrigido para uma propriedade específica, com seus próprios dados e validações.' },
  { term: 'Rotina', definition: 'Receita de cálculo pronta, construída e mantida no back-end. O painel só escolhe, aciona e acompanha; ela também roda pela linha de comando.' },
]

/** Sobre: explica o Delta Chemistry de forma geral, sem depender de dados ao vivo. */
export function AboutView() {
  return (
    <div className="view">
      <PageHero
        eyebrow="Sobre · Delta Chemistry" title="Métodos rápidos," accent="corrigidos por IA."
        lead="O Delta Chemistry é um hub de rotinas de cálculo de propriedades de moléculas. Ele usa métodos rápidos e corrige o erro deles com machine learning, aprendido a partir de dados de métodos mais precisos."
      />

      <section className="panel" aria-label="A ideia">
        <SectionHeading eyebrow="A ideia" title="Rapidez com precisão, em três passos">
          Em geral é preciso escolher entre um cálculo rápido e pouco preciso ou um cálculo preciso e demorado. O hub tenta reduzir essa diferença sem abandonar o método rápido.
        </SectionHeading>
        <div className="about-notes">
          <StoryNote number="1" title="Calcular rápido">O método rápido estima a propriedade da molécula em segundos, com algum erro.</StoryNote>
          <StoryNote number="2" title="Aprender o erro">Com moléculas calculadas também pelo método preciso, a IA aprende Δ = referência − rápido.</StoryNote>
          <StoryNote number="3" title="Corrigir">Para uma molécula nova, soma-se ao valor rápido a correção Δ prevista pelo modelo.</StoryNote>
        </div>
        <DeltaExplainer />
      </section>

      <section className="panel" aria-label="Como o hub se organiza">
        <SectionHeading eyebrow="Como o hub se organiza" title="Resumos e módulos" />
        <ul className="about-list">
          <li><strong>Resumos</strong><span>O compilado das propriedades obtidas pelas rotinas.</span></li>
          <li><strong>Módulos</strong><span>Cada propriedade tem o seu módulo na barra lateral, com dados, treino e validação próprios.</span></li>
        </ul>
        <h3 className="about-subtitle">Módulos disponíveis</h3>
        <ul className="about-list">
          {MODULES.map((module) => (
            <li key={module.key}><strong>{module.name}</strong><span>{module.property}. {module.description}</span></li>
          ))}
        </ul>
        <Callout label="Escopo atual" tone="warning">
          Hoje só existe o módulo de entalpia de formação no estado padrão. Outras propriedades só entram depois que houver dados, implementação e validação para elas.
        </Callout>
      </section>

      <Glossary terms={GLOSSARY} />
    </div>
  )
}

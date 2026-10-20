import { PageHero, SectionHeading } from './StoryVisuals'

const STEPS = [
  { title: 'Selecionar dados', text: 'Só entram moléculas com status OK e qualidade CBS OK; o filtro vem antes da divisão treino/teste.' },
  { title: 'Calcular o Δ', text: 'Δ = ΔHf CBS − ΔHf PM7. O modelo aprende o erro do PM7, não a energia inteira.' },
  { title: 'Ajustar o conjunto', text: 'KRR (similaridade suave) e XGBoost (interações não lineares), com média ponderada.' },
  { title: 'Avaliar no teste', text: 'Estatísticas calculadas só em moléculas que o modelo nunca viu, com o pré-processamento ajustado apenas no treino.' },
  { title: 'Portão de qualidade', text: 'O modelo só é aceito se MAE ≤ 3,5 kcal/mol, R² ≥ 0,97 e RMSE ≤ 5,0 kcal/mol.' },
]

export function TrainView() {
  return (
    <div className="view">
      <PageHero
        variant="control" eyebrow="Treinar modelo · Delta Learning" title="Ensine a IA a corrigir" accent="o cálculo rápido."
        lead="Aqui ficará o treino do modelo que soma ao PM7 a correção aprendida com os cálculos CBS. Por enquanto o treino roda pela linha de comando; esta tela explica o caminho."
      />
      <section className="panel" aria-label="Etapas do treino">
        <SectionHeading eyebrow="Como o treino funciona" title="Cinco etapas, sem vazamento de dados" badge={<span className="badge">Em breve no painel</span>} />
        <ol className="train-steps">
          {STEPS.map((step, index) => (
            <li key={step.title}><span aria-hidden="true">{index + 1}</span><div><strong>{step.title}</strong><p>{step.text}</p></div></li>
          ))}
        </ol>
        <p className="train-note">O painel aciona e acompanha as rotinas: quando o Delta Learning existir no back-end, o treino será configurado e acompanhado daqui, como em Calcular PM7.</p>
      </section>
    </div>
  )
}

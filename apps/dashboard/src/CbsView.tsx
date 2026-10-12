import type { DatasetProfile } from './api'
import { DataTable } from './DataTable'
import { findHistogram } from './helpers'
import { CompositionCard, EnthalpyCard, SizeCard } from './ProfileCharts'
import { ErrorNotice } from './Resource'
import {
  Callout, ChartCard, Explainer, Glossary, MetricCard, PageHero, ProportionBar, SectionHeading,
} from './StoryVisuals'
import { formatNumber, series } from './theme'
import { useApi } from './useApi'

export function MissingTable({ message, retry, file }: { message: string; retry: () => void; file: string }) {
  const missing = message.includes('file not found')
  return missing ? (
    <section className="panel missing-table" aria-label="Arquivo não encontrado">
      <span className="missing-table__icon" aria-hidden="true">🗂</span>
      <div>
        <h2>O arquivo <code>{file}</code> não está neste computador</h2>
        <p>Os CSVs grandes ficam fora do Git (estão no <code>.gitignore</code>). Copie o arquivo para a pasta <code>data/</code> do repositório e clique em atualizar.</p>
        <button type="button" className="ghost-button" onClick={retry}>Tentar novamente</button>
      </div>
    </section>
  ) : <ErrorNotice message={message} retry={retry} />
}

const CBS_GLOSSARY = [
  { term: 'CBS (Complete Basis Set)', definition: 'Família de métodos quânticos que extrapola o resultado para uma “base completa”, eliminando boa parte do erro numérico. Muito preciso e muito caro.' },
  { term: 'B3LYP', definition: 'Método DFT mais barato, usado como etapa intermediária. A diferença CBS − B3LYP mostra quanto a correção de alto nível mudou o valor.' },
  { term: 'Multiplicidade', definition: '1 = todos os elétrons emparelhados (molécula comum); 2 = um elétron desemparelhado (radical).' },
  { term: 'Átomos pesados', definition: 'Todos os átomos exceto hidrogênio. Medida simples do tamanho da molécula.' },
  { term: 'Conformação', definition: 'Uma das formas 3D que a mesma molécula pode assumir girando suas ligações.' },
]

export function CbsView() {
  const { state, reload } = useApi<DatasetProfile>('/api/tables/cbs/profile')
  const profile = state.status === 'ready' ? state.data : null
  const correction = profile ? findHistogram(profile, 'cbs_b3') : undefined
  const multiplicity = profile?.categories.multiplicity ?? {}
  const duplicates = profile?.categories.smiles_count ?? {}
  const repeated = Object.entries(duplicates).filter(([count]) => count !== '1' && count !== 'unset').reduce((sum, [, n]) => sum + n, 0)
  const citation = typeof profile?.manifest?.source_citation === 'string' ? profile.manifest.source_citation : null
  return (
    <div className="view">
      <PageHero
        variant="cbs" eyebrow="03 · Tabela CBS" title="A referência calculada" accent="com rigor máximo."
        lead="Dezenas de milhares de moléculas com entalpia calculada por um método quântico de alta precisão: é o “gabarito” que ensina o modelo a corrigir o PM7."
        onRefresh={reload}
      />
      <Explainer
        title="O que é a tabela CBS?"
        intro="É um banco público de propriedades termoquímicas calculadas com métodos CBS, filtrado para moléculas que só contêm carbono, hidrogênio, oxigênio e nitrogênio. Medir tudo isso em laboratório seria impossível; calcular com CBS é caro, mas viável uma única vez."
        analogy="O PM7 é um rascunho rápido; o CBS é a versão revisada por um especialista. O modelo aprende a transformar um no outro."
        facts={[
          { icon: '🏛', title: 'Origem', text: citation ?? 'ThermoCBS — Dobbelaere et al., J. Cheminf. (2024), licença CC-BY 4.0.' },
          { icon: '🎯', title: 'Para que serve', text: 'Alvo do Δ-learning: o modelo aprende Δ = CBS − PM7 para cada molécula.' },
          { icon: '🔧', title: 'Tratamento', text: 'Filtro CHON e correção de unidades (cerca de 1% dos valores vinha em cal/mol).' },
        ]}
      />
      {state.status === 'loading' && <p role="status" className="loading-state">Carregando tabela CBS…</p>}
      {state.status === 'error' && <MissingTable message={state.message} retry={reload} file="data/thermo_cbs_chon.csv" />}
      {profile && (
        <>
          <div className="metric-grid metric-grid--four">
            <MetricCard label="Linhas na tabela" value={profile.rows.toLocaleString('pt-BR')} note="Cada linha é uma geometria calculada" accent icon="🧮" />
            <MetricCard label="Moléculas distintas" value={(profile.unique_smiles ?? 0).toLocaleString('pt-BR')} note="SMILES únicos" icon="⌬" />
            <MetricCard label="Camada fechada" value={(multiplicity['1'] ?? 0).toLocaleString('pt-BR')} note={`${(multiplicity['2'] ?? 0).toLocaleString('pt-BR')} radicais (multiplicidade 2)`} icon="◉" />
            <MetricCard label="Correção CBS típica" value={formatNumber(correction?.median ?? null, 1)} unit="kcal/mol" note="Mediana de CBS − B3LYP" icon="Δ" hint="Quanto o método caro muda o resultado do método intermediário." />
          </div>

          <SectionHeading eyebrow="Distribuições" title="Como a referência CBS se distribui">
            Os mesmos três olhares das outras tabelas — tamanho, elementos e energia — mais a correção que o CBS aplica.
          </SectionHeading>
          <div className="chart-grid">
            <EnthalpyCard histogram={findHistogram(profile, 'H298_cbs')} title="Quão estáveis são as moléculas?" subtitle="Entalpia de formação CBS (kcal/mol)." />
            <SizeCard profile={profile} color={series[2]} />
          </div>
          <div className="chart-grid">
            <CompositionCard profile={profile} />
            <EnthalpyCard
              histogram={correction} title="Quanto o CBS corrige o B3LYP?" color={series[1]}
              subtitle="Diferença CBS − B3LYP por molécula (kcal/mol)."
              refs={correction?.median != null ? [{ value: 0, label: 'zero' }, { value: correction.median, label: 'mediana' }] : []}
              read={<>Se o B3LYP já fosse perfeito, tudo estaria em zero. A correção típica é de <strong>{formatNumber(correction?.median ?? null, 1)} kcal/mol</strong> — prova de que métodos baratos erram de forma sistemática, e é exatamente esse tipo de erro que o Δ-learning aprende.</>}
            />
          </div>
          <div className="chart-grid">
            <ChartCard title="Spin e carga" subtitle="Quase todas as moléculas são neutras e de camada fechada.">
              <ProportionBar label="Multiplicidade" segments={Object.entries(multiplicity).map(([key, count], index) => ({
                key, label: key === '1' ? 'Camada fechada (1)' : key === '2' ? 'Radical (2)' : `Multiplicidade ${key}`, count, color: series[index % series.length],
              }))} />
              <ProportionBar label="Carga" segments={Object.entries(profile.categories.charge ?? {}).map(([key, count], index) => ({
                key, label: key === '0' ? 'Neutra (0)' : `Carga ${key}`, count, color: index === 0 ? series[0] : 'var(--status-neutral)',
              }))} />
            </ChartCard>
            <ChartCard title="Moléculas com mais de uma geometria" subtitle="Algumas moléculas aparecem com conformações distintas.">
              {Object.keys(duplicates).length ? (
                <>
                  <p className="big-number"><strong>{repeated.toLocaleString('pt-BR')}</strong> linhas pertencem a SMILES repetidos</p>
                  <Callout label="Por que importa">Na divisão treino/teste, a mesma molécula não pode cair dos dois lados — senão o teste vira cola.</Callout>
                </>
              ) : <p className="empty-state">A coluna smiles_count não existe neste CSV.</p>}
            </ChartCard>
          </div>

          <section className="panel" aria-label="Tabela CBS completa">
            <SectionHeading eyebrow="Consulta" title="A tabela completa">Pesquise por ID ou SMILES; ordene por qualquer coluna.</SectionHeading>
            <DataTable table="cbs" title="Referência CBS" searchHint="Ex.: cbs_00012 ou CCO" defaultSort={{ column: 'mol_id', desc: false }} />
          </section>
        </>
      )}
      <Glossary terms={CBS_GLOSSARY} />
    </div>
  )
}

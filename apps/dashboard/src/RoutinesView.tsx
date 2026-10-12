import { useState } from 'react'
import type { RoutineCatalog, RoutineParameter, RoutineProgram } from './api'
import { ErrorNotice } from './Resource'
import {
  defaults, loadRoutines, newRoutine, parseMolecules, pipelineProblems, PRESETS, presetOf, routineProblems, saveRoutines,
} from './routines'
import type { ParamSet, ParamValue, RoutineSpec } from './routines'
import { Callout, PageHero } from './StoryVisuals'
import { useApi } from './useApi'

const STEPS = ['Nome', 'Propriedade', 'Pipeline', 'Configuração', 'Saída'] as const

const PROGRAM_ICONS: Record<string, string> = { rdkit: '⌬', crest: '⟳', confpass: '⇣', mopac: '⚡', smiles: '✎', name: '🔎', manual: '⌨', spreadsheet: '▦', value: '◎' }

function ParameterField({ parameter, value, onChange }: { parameter: RoutineParameter; value: ParamValue; onChange: (value: ParamValue) => void }) {
  const id = `param-${parameter.key}`
  const changed = value !== parameter.default
  const label = (
    <span className="param__label">
      {parameter.label}{parameter.unit && <small> ({parameter.unit})</small>}
      {changed && <button type="button" className="param__reset" onClick={() => onChange(parameter.default)} title="Voltar ao padrão">padrão: {String(parameter.default ?? 'vazio')}</button>}
    </span>
  )
  if (parameter.kind === 'bool') {
    return (
      <div className="param param--bool">
        <label className="switch"><input id={id} type="checkbox" checked={Boolean(value)} onChange={(event) => onChange(event.target.checked)} /><span aria-hidden="true" />{label}</label>
        {parameter.help && <p className="param__help">{parameter.help}</p>}
      </div>
    )
  }
  return (
    <div className="param">
      <label htmlFor={id}>{label}</label>
      {parameter.kind === 'choice' ? (
        <select id={id} value={String(value ?? '')} onChange={(event) => onChange(event.target.value)}>
          {parameter.choices.map((choice) => <option key={choice} value={choice}>{choice}</option>)}
        </select>
      ) : (
        <input
          id={id} type={parameter.kind === 'text' ? 'text' : 'number'} value={value === null ? '' : String(value)}
          step={parameter.kind === 'int' ? 1 : 'any'} min={parameter.minimum ?? undefined} max={parameter.maximum ?? undefined}
          onChange={(event) => {
            const raw = event.target.value
            if (parameter.kind === 'text') onChange(raw)
            else onChange(raw === '' ? null : parameter.kind === 'int' ? Math.trunc(Number(raw)) : Number(raw))
          }}
        />
      )}
      {parameter.help && <p className="param__help">{parameter.help}</p>}
    </div>
  )
}

function ProgramForm({ program, values, onChange }: { program: RoutineProgram; values: ParamSet; onChange: (values: ParamSet) => void }) {
  return (
    <fieldset className="program-form">
      <legend>
        <span className="program-form__icon" aria-hidden="true">{PROGRAM_ICONS[program.key] ?? '•'}</span>
        {program.label}
        {program.experimental && <span className="pill pill--warning">experimental</span>}
        {program.available === false && <span className="pill pill--critical">não encontrado no PATH</span>}
        {program.available === true && <span className="pill pill--ok">instalado</span>}
      </legend>
      <p className="program-form__text">{program.description}</p>
      {program.parameters.length === 0 ? <p className="subtle">Sem parâmetros a ajustar.</p> : (
        <div className="param-grid">
          {program.parameters.map((parameter) => (
            <ParameterField key={parameter.key} parameter={parameter} value={values[parameter.key] ?? parameter.default}
              onChange={(value) => onChange({ ...values, [parameter.key]: value })} />
          ))}
        </div>
      )}
      <button type="button" className="ghost-button ghost-button--small" onClick={() => onChange(defaults(program))}>Restaurar padrões</button>
    </fieldset>
  )
}

function Wizard({ catalog, initial, onSave, onCancel }: {
  catalog: RoutineCatalog
  initial: RoutineSpec
  onSave: (spec: RoutineSpec) => void
  onCancel: () => void
}) {
  const [spec, setSpec] = useState(initial)
  const [step, setStep] = useState(1)
  const [moleculeText, setMoleculeText] = useState(initial.input.molecules.join('\n'))
  const update = (patch: Partial<RoutineSpec>) => setSpec((current) => ({ ...current, ...patch }))
  const problems = routineProblems(spec)
  const programs = Object.fromEntries(catalog.programs.map((program) => [program.key, program]))
  const inputMode = catalog.input_modes.find((mode) => mode.key === spec.input.mode)
  const preset = presetOf(spec.steps)
  const pipeIssues = pipelineProblems(spec.steps)

  const move = (index: number, offset: number) => {
    const next = [...spec.steps]
    const target = index + offset
    if (target < 0 || target >= next.length) return
    ;[next[index], next[target]] = [next[target], next[index]]
    update({ steps: next })
  }
  const toggle = (key: string) => {
    if (spec.steps.includes(key)) update({ steps: spec.steps.filter((item) => item !== key) })
    else {
      const order = ['rdkit', 'crest', 'confpass', 'mopac']
      update({ steps: [...spec.steps, key].sort((a, b) => order.indexOf(a) - order.indexOf(b)) })
    }
  }
  const download = () => {
    const blob = new Blob([JSON.stringify(spec, null, 2)], { type: 'application/json' })
    const url = URL.createObjectURL(blob)
    const anchor = document.createElement('a')
    anchor.href = url
    anchor.download = `${spec.name.trim().replace(/[^\w-]+/g, '_') || 'rotina'}.json`
    anchor.click()
    URL.revokeObjectURL(url)
  }

  return (
    <div className="wizard">
      <ol className="stepper" aria-label="Passos da rotina">
        {STEPS.map((label, index) => {
          const number = index + 1
          const hasProblem = problems.some((problem) => problem.step === number)
          return (
            <li key={label}>
              <button type="button" aria-current={step === number ? 'step' : undefined} className={`stepper__item${number < step ? ' stepper__item--done' : ''}${hasProblem && number < step ? ' stepper__item--warn' : ''}`} onClick={() => setStep(number)}>
                <span className="stepper__number">{number < step && !hasProblem ? '✓' : number}</span>
                <span>{label}</span>
              </button>
            </li>
          )
        })}
        <li>
          <button type="button" aria-current={step === 6 ? 'step' : undefined} className="stepper__item" onClick={() => setStep(6)}>
            <span className="stepper__number">★</span><span>Revisão</span>
          </button>
        </li>
      </ol>

      <div className="wizard__body">
        {step === 1 && (
          <section aria-label="Passo 1: nome">
            <h3>Como vamos chamar esta rotina?</h3>
            <p className="wizard__lead">Um nome curto ajuda a reencontrar a rotina no menu.</p>
            <div className="param"><label htmlFor="routine-name">Nome da rotina</label>
              <input id="routine-name" value={spec.name} onChange={(event) => update({ name: event.target.value })} placeholder="Ex.: Álcoois C2–C6 · PM7 completo" /></div>
            <div className="param"><label htmlFor="routine-description">Descrição (opcional)</label>
              <textarea id="routine-description" rows={3} value={spec.description} onChange={(event) => update({ description: event.target.value })} placeholder="Para que serve esta rotina?" /></div>
          </section>
        )}

        {step === 2 && (
          <section aria-label="Passo 2: propriedade">
            <h3>Qual propriedade calcular?</h3>
            <p className="wizard__lead">As propriedades disponíveis vêm do back-end. Hoje o sistema calcula uma.</p>
            <div className="choice-grid">
              {catalog.properties.map((property) => (
                <button key={property.key} type="button" disabled={!property.available} aria-pressed={spec.property === property.key}
                  className="choice-card" onClick={() => update({ property: property.key })}>
                  <span className="choice-card__icon" aria-hidden="true">ΔH</span>
                  <strong>{property.label}</strong>
                  <span>{property.description}</span>
                  <small>Unidade: {property.unit}</small>
                </button>
              ))}
              <div className="choice-card choice-card--soon" aria-hidden="true"><span className="choice-card__icon">+</span><strong>Novas propriedades</strong><span>Entropia, capacidade calorífica… quando o pipeline suportar.</span></div>
            </div>
          </section>
        )}

        {step === 3 && (
          <section aria-label="Passo 3: pipeline">
            <h3>3.1 · Como as moléculas entram?</h3>
            <div className="choice-grid choice-grid--three">
              {catalog.input_modes.map((mode) => (
                <button key={mode.key} type="button" aria-pressed={spec.input.mode === mode.key} className="choice-card"
                  onClick={() => update({ input: { ...spec.input, mode: mode.key, params: defaults(mode) } })}>
                  <span className="choice-card__icon" aria-hidden="true">{PROGRAM_ICONS[mode.key] ?? '•'}</span>
                  <strong>{mode.label}</strong>
                  <span>{mode.description}</span>
                </button>
              ))}
            </div>
            <div className="param">
              <label htmlFor="routine-molecules">{spec.input.mode === 'name' ? 'Nomes das moléculas (um por linha)' : 'Estruturas (uma por linha)'}</label>
              <textarea id="routine-molecules" rows={5} value={moleculeText} className="mono"
                placeholder={spec.input.mode === 'name' ? 'ethanol\nacetone\nbenzene' : 'CCO\nCC(C)=O\nc1ccccc1'}
                onChange={(event) => { setMoleculeText(event.target.value); update({ input: { ...spec.input, molecules: parseMolecules(event.target.value) } }) }} />
              <p className="param__help">{spec.input.molecules.length.toLocaleString('pt-BR')} molécula(s) informada(s).</p>
            </div>

            <h3>3.2 · Quais programas rodam, e em que ordem?</h3>
            <div className="choice-grid choice-grid--three">
              {Object.entries(PRESETS).map(([key, item]) => (
                <button key={key} type="button" aria-pressed={preset === key} className="choice-card" onClick={() => update({ steps: [...item.steps] })}>
                  <strong>{item.label}</strong><span>{item.text}</span>
                  <small className="mono">{item.steps.join(' → ')}</small>
                </button>
              ))}
              <div className={`choice-card choice-card--static${preset === 'custom' ? ' choice-card--active' : ''}`}><strong>Personalizado</strong><span>Ligue, desligue e reordene os programas abaixo.</span></div>
            </div>
            <div className="flow-builder">
              <div className="flow-strip" aria-label="Ordem atual do pipeline">
                <span className="flow-chip flow-chip--input">{PROGRAM_ICONS[spec.input.mode]} {inputMode?.label ?? spec.input.mode}</span>
                {spec.steps.map((key, index) => (
                  <span key={key} className="flow-chip-wrap">
                    <span className="flow-arrow" aria-hidden="true">→</span>
                    <span className="flow-chip">
                      {PROGRAM_ICONS[key]} {programs[key]?.label ?? key}
                      <span className="flow-chip__actions">
                        <button type="button" onClick={() => move(index, -1)} disabled={index === 0} aria-label={`Mover ${key} para antes`}>‹</button>
                        <button type="button" onClick={() => move(index, 1)} disabled={index === spec.steps.length - 1} aria-label={`Mover ${key} para depois`}>›</button>
                      </span>
                    </span>
                  </span>
                ))}
                <span className="flow-arrow" aria-hidden="true">→</span>
                <span className="flow-chip flow-chip--output">ΔHf°</span>
              </div>
              <div className="toggle-row">
                {catalog.programs.map((program) => (
                  <label key={program.key} className="switch">
                    <input type="checkbox" checked={spec.steps.includes(program.key)} onChange={() => toggle(program.key)} />
                    <span aria-hidden="true" />{program.label}
                  </label>
                ))}
              </div>
              {pipeIssues.length > 0
                ? <ul className="issue-list" role="alert">{pipeIssues.map((issue) => <li key={issue}>{issue}</li>)}</ul>
                : <p className="ok-line">✓ Ordem válida para o pipeline Semi-Imperium.</p>}
            </div>
          </section>
        )}

        {step === 4 && (
          <section aria-label="Passo 4: configuração">
            <h3>Ajuste cada programa</h3>
            <p className="wizard__lead">Os valores já vêm com os padrões do Semi-Imperium. Mude só o que precisar — cada alteração mostra o padrão para voltar.</p>
            {inputMode && inputMode.parameters.length > 0 && (
              <ProgramForm program={inputMode} values={spec.input.params} onChange={(params) => update({ input: { ...spec.input, params } })} />
            )}
            {spec.steps.map((key, index) => programs[key] && (
              <div key={key} className="program-step"><span className="program-step__index">4.{index + 1}</span>
                <ProgramForm program={programs[key]} values={spec.programs[key] ?? {}} onChange={(params) => update({ programs: { ...spec.programs, [key]: params } })} />
              </div>
            ))}
          </section>
        )}

        {step === 5 && (
          <section aria-label="Passo 5: saída">
            <h3>Como você quer o resultado?</h3>
            <div className="choice-grid choice-grid--three">
              {catalog.outputs.map((output) => {
                const active = spec.outputs.kinds.includes(output.key)
                return (
                  <button key={output.key} type="button" aria-pressed={active} className="choice-card"
                    onClick={() => update({ outputs: { ...spec.outputs, kinds: active ? spec.outputs.kinds.filter((k) => k !== output.key) : [...spec.outputs.kinds, output.key] } })}>
                    <span className="choice-card__icon" aria-hidden="true">{PROGRAM_ICONS[output.key] ?? '•'}</span>
                    <strong>{output.label}</strong><span>{output.description}</span>
                  </button>
                )
              })}
            </div>
            {catalog.outputs.filter((output) => spec.outputs.kinds.includes(output.key) && output.parameters.length > 0).map((output) => (
              <ProgramForm key={output.key} program={output} values={spec.outputs.params[output.key] ?? {}}
                onChange={(params) => update({ outputs: { ...spec.outputs, params: { ...spec.outputs.params, [output.key]: params } } })} />
            ))}
          </section>
        )}

        {step === 6 && (
          <section aria-label="Revisão">
            <h3>Revisão da rotina</h3>
            {problems.length > 0 ? (
              <ul className="issue-list" role="alert">{problems.map((problem) => (
                <li key={problem.text}><button type="button" className="link-button" onClick={() => setStep(problem.step)}>Passo {problem.step}</button>: {problem.text}</li>
              ))}</ul>
            ) : <p className="ok-line">✓ Tudo pronto. A rotina pode ser salva e exportada.</p>}
            <dl className="review-list">
              <div><dt>Nome</dt><dd>{spec.name || '—'}</dd></div>
              <div><dt>Propriedade</dt><dd>{catalog.properties.find((item) => item.key === spec.property)?.label ?? '—'}</dd></div>
              <div><dt>Entrada</dt><dd>{inputMode?.label} · {spec.input.molecules.length} molécula(s)</dd></div>
              <div><dt>Pipeline</dt><dd className="mono">{spec.steps.join(' → ')}</dd></div>
              <div><dt>Saída</dt><dd>{spec.outputs.kinds.map((kind) => catalog.outputs.find((o) => o.key === kind)?.label ?? kind).join(' + ') || '—'}</dd></div>
            </dl>
            <details className="json-preview"><summary>Ver especificação (JSON enviado ao back-end)</summary><pre>{JSON.stringify(spec, null, 2)}</pre></details>
            <Callout label="Disparo" tone="warning">
              O botão “Disparar” será ativado quando o back-end de jobs (fila, status e log) estiver pronto. Até lá, salve a rotina ou exporte o JSON.
            </Callout>
          </section>
        )}
      </div>

      <div className="wizard__footer">
        <button type="button" className="ghost-button" onClick={onCancel}>Cancelar</button>
        <div className="wizard__nav">
          {step > 1 && <button type="button" className="ghost-button" onClick={() => setStep(step - 1)}>← Voltar</button>}
          {step < 6 && <button type="button" className="primary-button" onClick={() => setStep(step + 1)}>Próximo →</button>}
          {step === 6 && (
            <>
              <button type="button" className="ghost-button" onClick={download}>Exportar JSON</button>
              <button type="button" className="primary-button" disabled={problems.length > 0} onClick={() => onSave({ ...spec, updated_at: new Date().toISOString() })}>Salvar rotina</button>
              <button type="button" className="primary-button primary-button--launch" disabled={!catalog.dispatch_enabled} title="Disponível quando o back-end de jobs estiver pronto">▶ Disparar</button>
            </>
          )}
        </div>
      </div>
    </div>
  )
}

export function RoutinesView() {
  const { state, reload } = useApi<RoutineCatalog>('/api/routines/catalog')
  const [routines, setRoutines] = useState<RoutineSpec[]>(() => loadRoutines())
  const [editing, setEditing] = useState<RoutineSpec | null>(null)
  const [notice, setNotice] = useState<string | null>(null)
  const catalog = state.status === 'ready' ? state.data : null

  const persist = (next: RoutineSpec[], message: string) => {
    setRoutines(next)
    setNotice(saveRoutines(next) ? message : `${message} Atenção: o navegador bloqueou o armazenamento local; a rotina some ao recarregar.`)
  }

  return (
    <div className="view">
      <PageHero
        variant="control" eyebrow="Calcular · Painel de controle" title="Monte a rotina," accent="o back-end executa."
        lead="O dashboard não calcula nada: ele descreve o que deve ser feito e dispara o processo. Monte rotinas em cinco passos com os programas e parâmetros que o Semi-Imperium aceita."
        art={<div className="control-art" aria-hidden="true"><span>⌬</span><i /><span>⟳</span><i /><span>⇣</span><i /><span>⚡</span></div>}
      />
      {state.status === 'loading' && <p role="status" className="loading-state">Carregando programas disponíveis…</p>}
      {state.status === 'error' && <ErrorNotice message={state.message} retry={reload} />}
      {notice && <p role="status" className="toast">{notice} <button type="button" onClick={() => setNotice(null)} aria-label="Fechar aviso">×</button></p>}
      {catalog && (
        <div className="routines-layout">
          <aside className="routine-menu" aria-label="Menu de rotinas">
            <div className="routine-menu__head"><h2>Rotinas</h2><span>{routines.length}</span></div>
            <button type="button" className="primary-button routine-menu__new" onClick={() => setEditing(newRoutine(catalog))}>＋ Criar nova rotina</button>
            {routines.length === 0 ? <p className="subtle">Nenhuma rotina salva ainda.</p> : (
              <ul>
                {routines.map((routine) => (
                  <li key={routine.id} className={editing?.id === routine.id ? 'routine-menu__item routine-menu__item--active' : 'routine-menu__item'}>
                    <button type="button" className="routine-menu__open" onClick={() => setEditing(routine)}>
                      <strong>{routine.name}</strong>
                      <span>{routine.input.molecules.length} molécula(s) · {routine.steps.join(' → ')}</span>
                    </button>
                    <div className="routine-menu__actions">
                      <button type="button" onClick={() => setEditing({ ...routine, id: `routine-${Date.now().toString(36)}`, name: `${routine.name} (cópia)` })}>Duplicar</button>
                      <button type="button" onClick={() => { persist(routines.filter((item) => item.id !== routine.id), `Rotina “${routine.name}” removida.`); if (editing?.id === routine.id) setEditing(null) }}>Excluir</button>
                    </div>
                  </li>
                ))}
              </ul>
            )}
            <div className="tool-status">
              <h3>Programas nesta máquina</h3>
              <ul>
                {catalog.programs.filter((program) => program.executable).map((program) => (
                  <li key={program.key}><span className={`dot dot--${program.available ? 'ok' : 'missing'}`} aria-hidden="true" />{program.label}<small>{program.available ? 'instalado' : 'não encontrado'}</small></li>
                ))}
              </ul>
            </div>
          </aside>
          <section className="panel routine-main" aria-label="Editor de rotina">
            {editing ? (
              <Wizard key={editing.id} catalog={catalog} initial={editing} onCancel={() => setEditing(null)}
                onSave={(spec) => {
                  const exists = routines.some((item) => item.id === spec.id)
                  persist(exists ? routines.map((item) => item.id === spec.id ? spec : item) : [...routines, spec], `Rotina “${spec.name}” salva.`)
                  setEditing(null)
                }} />
            ) : (
              <div className="routine-empty">
                <span aria-hidden="true">🧭</span>
                <h2>Escolha uma rotina ou crie uma nova</h2>
                <ol>
                  <li><strong>Nome</strong> — como a rotina aparece no menu.</li>
                  <li><strong>Propriedade</strong> — o que calcular (hoje: ΔHf°).</li>
                  <li><strong>Pipeline</strong> — como as moléculas entram e quais programas rodam.</li>
                  <li><strong>Configuração</strong> — parâmetros de cada programa.</li>
                  <li><strong>Saída</strong> — planilha, valor na tela ou ambos.</li>
                </ol>
                <button type="button" className="primary-button" onClick={() => setEditing(newRoutine(catalog))}>＋ Criar nova rotina</button>
              </div>
            )}
          </section>
        </div>
      )}
    </div>
  )
}

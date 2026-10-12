import { useEffect, useState } from 'react'
import type { Cell, TableKey, TablePage } from './api'
import { COLUMN_INFO } from './helpers'
import { ErrorNotice } from './Resource'
import { formatNumber, formatSigned } from './theme'
import { useApi } from './useApi'

const PAGE_SIZE = 25

function render(column: string, value: Cell) {
  if (value === null) return <span className="cell-empty">—</span>
  if (typeof value === 'boolean') return value ? <span className="pill pill--flag">sim</span> : <span className="cell-empty">não</span>
  if (typeof value === 'number') return COLUMN_INFO[column]?.signed ? formatSigned(value) : formatNumber(value)
  if (column === 'status' || column.endsWith('_status') || column === 'cbs_quality_flag') {
    return <span className={`pill pill--${value.toLowerCase()}`}>{value}</span>
  }
  return value
}

export function DataTable({ table, title, searchHint, defaultSort }: {
  table: TableKey
  title: string
  searchHint: string
  defaultSort?: { column: string; desc: boolean }
}) {
  const [draft, setDraft] = useState('')
  const [query, setQuery] = useState('')
  const [page, setPage] = useState(0)
  const [sort, setSort] = useState(defaultSort ?? null)

  useEffect(() => {
    const timer = setTimeout(() => { setQuery(draft.trim()); setPage(0) }, 250)
    return () => clearTimeout(timer)
  }, [draft])

  const params = new URLSearchParams({ offset: String(page * PAGE_SIZE), limit: String(PAGE_SIZE) })
  if (query) params.set('q', query)
  if (sort) { params.set('sort', sort.column); params.set('desc', String(sort.desc)) }
  const { state, reload } = useApi<TablePage>(`/api/tables/${table}/rows?${params.toString()}`)

  const toggleSort = (column: string) => {
    setSort((current) => current?.column === column ? { column, desc: !current.desc } : { column, desc: false })
    setPage(0)
  }

  return (
    <section className="data-table" aria-label={title}>
      <div className="data-table__toolbar">
        <div>
          <h3>{title}</h3>
          {state.status === 'ready' && (
            <p>{state.data.matched === state.data.total
              ? `${state.data.total.toLocaleString('pt-BR')} linhas`
              : `${state.data.matched.toLocaleString('pt-BR')} de ${state.data.total.toLocaleString('pt-BR')} linhas correspondem à busca`}
            {' · '}clique no cabeçalho para ordenar</p>
          )}
        </div>
        <label className="search-field">
          <span className="sr-only">Buscar em {title}</span>
          <span aria-hidden="true">⌕</span>
          <input type="search" value={draft} onChange={(event) => setDraft(event.target.value)} placeholder={searchHint} />
        </label>
      </div>
      {state.status === 'loading' && <p role="status" className="loading-state">Carregando linhas…</p>}
      {state.status === 'error' && <ErrorNotice message={state.message} retry={reload} />}
      {state.status === 'ready' && (state.data.rows.length === 0 ? <p className="empty-state">Nenhuma linha corresponde à busca.</p> : (
        <>
          <div className="table-wrap">
            <table>
              <caption className="sr-only">{title}</caption>
              <thead>
                <tr>{state.data.columns.map((column) => {
                  const active = sort?.column === column
                  return (
                    <th key={column} scope="col" aria-sort={active ? (sort.desc ? 'descending' : 'ascending') : 'none'} title={COLUMN_INFO[column]?.help}>
                      <button type="button" className="sort-button" onClick={() => toggleSort(column)}>
                        {COLUMN_INFO[column]?.label ?? column}
                        <span aria-hidden="true">{active ? (sort.desc ? '↓' : '↑') : '↕'}</span>
                      </button>
                    </th>
                  )
                })}</tr>
              </thead>
              <tbody>{state.data.rows.map((row, index) => (
                <tr key={`${String(row[state.data.columns[0]])}-${index}`}>
                  {state.data.columns.map((column) => (
                    <td key={column} className={typeof row[column] === 'number' ? 'num-cell' : column === 'smiles' ? 'smiles' : column.endsWith('_id') ? 'id-cell' : undefined}>
                      {render(column, row[column])}
                    </td>
                  ))}
                </tr>
              ))}</tbody>
            </table>
          </div>
          <div className="pagination">
            <span>Mostrando {state.data.offset + 1}–{state.data.offset + state.data.rows.length} de {state.data.matched.toLocaleString('pt-BR')}</span>
            <button type="button" disabled={page === 0} onClick={() => setPage((value) => value - 1)}>Anterior</button>
            <button type="button" disabled={state.data.offset + state.data.rows.length >= state.data.matched} onClick={() => setPage((value) => value + 1)}>Próxima</button>
          </div>
        </>
      ))}
    </section>
  )
}

import { render, screen } from '@testing-library/react'
import { afterEach, expect, test } from 'vitest'
import { cleanup } from '@testing-library/react'
import { App } from './App'

afterEach(cleanup)

test('monta o scaffold do dashboard', () => {
  render(<App />)
  expect(screen.getByRole('heading', { name: 'Grimperium' })).toBeInTheDocument()
})

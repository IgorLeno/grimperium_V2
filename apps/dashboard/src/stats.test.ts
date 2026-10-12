import { expect, test } from 'vitest'
import type { ComparisonPoint } from './api'
import { niceTicks } from './helpers'
import { pipelineProblems, presetOf } from './routines'
import { byHeavyAtoms, bySource, cumulativeAccuracy, histogram, summarize, thresholdFor, tierCounts } from './stats'

function point(error: number, nheavy = 3, source = 'atct'): ComparisonPoint {
  return { exp_id: `e${error}`, mol_id: 'm', smiles: 'C', nheavy, source, H298_exp: 0, uncertainty: null,
    H298_cbs: error, cbs_minus_exp: error, conflict_flag: false, low_priority_only: false, validation_only: false }
}

test('summarize computes MAE, bias, RMSE and chemical-accuracy share', () => {
  const stats = summarize([point(1), point(-3), point(0.5)])
  expect(stats.mae).toBeCloseTo(1.5)
  expect(stats.bias).toBeCloseTo(-0.5)
  expect(stats.rmse).toBeCloseTo(Math.sqrt((1 + 9 + 0.25) / 3))
  expect(stats.withinChemical).toBeCloseTo(2 / 3)
  expect(summarize([]).mae).toBeNull()
})

test('tiers use inclusive upper bounds', () => {
  const counts = tierCounts([point(1), point(-5), point(10), point(10.01)])
  expect(counts.map((item) => item.count)).toEqual([1, 1, 1, 1])
})

test('histogram counts every value either in a bin or in the clipped tails', () => {
  const values = [...Array.from({ length: 200 }, (_, i) => (i % 10) - 5), 400]
  const result = histogram(values, 1, 0.005)
  const binned = result.bins.reduce((sum, bin) => sum + bin.count, 0)
  expect(binned + result.below + result.above).toBe(values.length)
  expect(result.above).toBe(1)
  expect(Number.isInteger(result.bins[0].start)).toBe(true)
})

test('cumulative accuracy reaches 100% and thresholdFor matches it', () => {
  const points = [point(0.2), point(-1.5), point(4)]
  const curve = cumulativeAccuracy(points, 0.5)
  expect(curve[curve.length - 1].share).toBe(100)
  expect(curve.find((item) => item.threshold === 1)?.share).toBeCloseTo(100 / 3)
  expect(thresholdFor(points, 66)).toBeCloseTo(1.5)
})

test('grouping by source and by size', () => {
  const points = [point(1, 2, 'a'), point(-3, 2, 'a'), point(2, 5, 'b'), ...Array.from({ length: 5 }, () => point(1, 5, 'c'))]
  expect(bySource(points, 5).map((item) => item.key)).toEqual(['c'])
  expect(byHeavyAtoms(points).map((item) => [item.nheavy, item.n])).toEqual([[2, 2], [5, 6]])
})

test('niceTicks picks round steps', () => {
  expect(niceTicks(-200, 160, 6)).toEqual([-200, -150, -100, -50, 0, 50, 100, 150])
  expect(niceTicks(0, 1, 5)).toEqual([0, 0.2, 0.4, 0.6, 0.8, 1])
})

test('pipeline rules follow the Semi-Imperium order', () => {
  expect(pipelineProblems(['rdkit', 'crest', 'confpass', 'mopac'])).toEqual([])
  expect(pipelineProblems(['rdkit', 'mopac'])).toEqual([])
  expect(pipelineProblems(['rdkit', 'confpass', 'mopac'])).toContain('A seleção de conformações precisa do CREST antes dela.')
  expect(pipelineProblems(['crest', 'rdkit', 'mopac'])).toContain('O RDKit deve ser o primeiro passo.')
  expect(presetOf(['rdkit', 'mopac'])).toBe('direct_pm7')
  expect(presetOf(['rdkit', 'crest', 'mopac'])).toBe('custom')
})

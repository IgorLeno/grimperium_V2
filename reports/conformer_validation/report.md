# Validação da seleção de conformeros

Gerado por `scripts/conformer_validation summarize`. Delta = `H298_cbs − ΔHf_PM7`, só mínimos verificados; DP amostral; MAD = desvio absoluto mediano. Conjuntos "comum": moléculas com valor em todos os braços.

| Braço | Moléculas | verified | saddle | unverified | failed | erro | n delta (comum) | DP delta (comum) | MAD delta (comum) | DP delta (todas) | Recuperação < 0,5 Å | RMSD final (mediana) | |X∩A| médio | Otimizações | FORCE | Custo s/molécula |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 200 | 199 | 0 | 0 | 0 | 1 | 195 | 5.984 | 3.760 | 6.028 | 0.800 | 0.466 | 6.945 | 1383 | 1383 | 56.558 |
| B | 200 | 195 | 0 | 0 | 0 | 5 | 195 | 5.985 | 3.760 | 5.985 | 0.815 | 0.473 | 6.087 | 1367 | 1367 | 55.501 |
| C | 200 | 195 | 0 | 0 | 0 | 5 | 195 | 5.985 | 3.760 | 5.985 | 0.815 | 0.471 | 6.056 | 1366 | 1366 | 56.292 |

## Decisão

**B** — Piso de recuperação 0.750 (A 0.800 - 0.05); elegíveis ['A', 'B', 'C']. Menor DP do delta 5.984 kcal/mol; dentro de 0.1 kcal/mol: ['A', 'B', 'C']. Empate: vence o mais barato, B.

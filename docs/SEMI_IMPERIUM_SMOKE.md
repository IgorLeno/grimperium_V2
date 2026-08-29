# Semi-Imperium: smoke científico real

Este registro distingue uma execução com binários reais dos testes normais com
doubles. Ele documenta o ambiente e o resultado observado; não é uma referência
termoquímica nem validação de exatidão.

## Execução de 2026-08-29

Comando:

```bash
SEMI_IMPERIUM_REAL_SMOKE=1 poetry run pytest \
  tests/integration/test_semi_imperium_real_smoke.py -v -s
```

Ambiente e configuração:

- molécula: água, SMILES `O`, carga 0, multiplicidade 1;
- CREST 3.0.2, GFN2, modo quick, uma thread;
- seleção `crest_energy_top_n`, Top-1;
- MOPAC 23.2.2, Hamiltoniano PM7;
- política `require_minimum`, sem tentativa de recuperação por deslocamento;
- Python 3.14.7, Linux;
- artefatos isolados no diretório temporário do pytest.

Resultado observado:

- teste: **PASS** em 7,54 s;
- ensemble CREST: 1 conformero;
- índice selecionado: 0;
- estado MOPAC: `minimum_verified`;
- `ΔHf°` provisória: -57,79986 kcal/mol;
- `ΔHf°` verificada: -57,79986 kcal/mol.

Esse PASS comprova somente que, neste ambiente, os adapters e os dois
executáveis completaram o pequeno percurso CREST → Energy Top-N → otimização
PM7 → `FORCE`, produzindo artefatos que os parsers do Semi-Imperium aceitaram.
Um único caso não mede erro do método, não valida CONFPASS, AM1 ou PM3, não
estabelece reprodutibilidade entre plataformas e não prova mínimo global.

Em ambientes sem CREST ou MOPAC, ou sem
`SEMI_IMPERIUM_REAL_SMOKE=1`, o teste é explicitamente `SKIPPED`. Nesse caso os
testes com doubles continuam úteis para contratos de software, mas não contam
como esta evidência científica real.

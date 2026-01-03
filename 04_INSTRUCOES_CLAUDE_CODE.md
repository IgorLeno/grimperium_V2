📝 INSTRUÇÕES: Como Usar o Prompt com Claude Code

🔄 Passo 1: Preparar o Prompt

Opção A: Copiar direto do arquivo
1. Abra no seu editor: 02_PROMPT_001_SCAFFOLDING_INICIAL.md
2. Select All (Ctrl+A ou Cmd+A)
3. Copy (Ctrl+C ou Cmd+C)

Opção B: Abrir e ler inline
Se preferir ler o prompt antes de colar, você tem o arquivo completo.

---

📌 Passo 2: Abrir Claude Code no Cursor

1. Abra o Cursor no seu projeto grimperium
2. Pressione Cmd+K (Mac) ou Ctrl+K (Windows/Linux)
   - Abre a "Command Palette" do Claude Code
3. Digite: @plan
   - Ou simplesmente comece digitando seu prompt

Alternativa: Chat Direto
1. Pressione Cmd+Shift+0 (Mac) ou Ctrl+Shift+0 (Windows)
2. Abre o painel de Chat do Claude

---

📌 Passo 3: Colar o Prompt

Opção A: Via Command Palette
1. Cmd+K → digitou @plan → Enter
2. Abre um painel de planejamento
3. Cole o prompt aqui (Cmd+V ou Ctrl+V)
4. Enter ou Send

Opção B: Via Chat Panel
1. Abrir chat (Cmd+Shift+0)
2. Cole o prompt inteiro na caixa de texto
3. Send

---

⏱️ Passo 4: Deixar o Claude Code Planejar

Claude Code vai automaticamente:

1. Ler o prompt
2. Executar @plan (se você usou @plan no início)
3. Quebrar em batches:
   - Batch 1: Estrutura de pastas
   - Batch 2: Configuração (Poetry, tox, pre-commit)
   - Batch 3: Módulos base
   - Batch 4: Testes e fixtures
   - Batch 5: Documentação

4. Você verá logs de progresso

5. Claude Code começa a executar batch por batch

---

👀 Passo 5: Acompanhar a Execução

Batch 1 (Estrutura) ~3 min
- Cria pastas: src/, tests/, docs/, .github/
- Cria arquivos vazios com comentários

Batch 2 (Configuração) ~5 min
- Cria pyproject.toml (Poetry)
- Cria tox.ini (multi-Python testing)
- Cria .pre-commit-config.yaml
- Cria .github/workflows/ci.yml
- Cria .gitignore, LICENSE

Batch 3 (Módulos) ~8 min
- Cria todos os .py em src/grimperium/
- Cada arquivo com imports + docstrings + stubs
- Nenhuma lógica real ainda

Batch 4 (Testes) ~5 min
- Cria tests/fixtures/mock_data.py com pytest fixtures
- Cria stubs em tests/unit/ e tests/integration/

Batch 5 (Docs) ~4 min
- Cria README.md com ASCII architecture
- Cria docs/*.md (architecture, guides)
- Cria CHANGELOG.md

Total esperado: ~25 minutos

✅ Sinais de Sucesso

Você verá logs como:
```
✓ Created: src/grimperium/__init__.py
✓ Created: pyproject.toml
✓ Created: tests/fixtures/mock_data.py
...
```

⚠️ Se houver erro

Claude Code vai reportar o erro. Opções:
1. @debug - para Claude Code debugar
2. @refactor - para refatorar abordagem
3. Aperte Esc para parar e pergunte manualmente

---

🧪 Passo 6: Validar o Resultado

Após Claude Code terminar, execute esses comandos no terminal:

1️⃣ Verificar estrutura

```bash
# Listar estrutura criada
ls -R src/grimperium/
ls -R tests/
ls -R docs/

# Deve conter ~25 arquivos
find . -type f -name "*.py" | wc -l
# Esperado: ~20-25
```

2️⃣ Instalar dependências

```bash
# Poetry install
poetry install

# Ativa environment
poetry shell
```

Esperado: Sem erros, todas as deps instaladas

3️⃣ Verificar imports

```bash
# Test se imports funcionam
python -c "from grimperium.models import BaseModel; print('✅ Imports OK')"
python -c "from grimperium.data import loader; print('✅ Data OK')"
python -c "from grimperium import api; print('✅ API OK')"
```

4️⃣ Executar linting

```bash
# Ruff
ruff check .
# Esperado: Maybe some warnings ok, no critical errors

# Black
black --check .
# Esperado: Sem reformatações necessárias

# Mypy
mypy src/
# Esperado: Sem erros críticos (warnings ok)
```

5️⃣ Executar testes

```bash
# Rodar testes
pytest tests/ -v

# Com coverage
pytest --cov=src/grimperium tests/
# Esperado: 100% pass (stubs são ok, alguns skip ok)
```

6️⃣ Tox (multi-Python)

```bash
# Se tiver 3.9, 3.10, 3.11, 3.12 instalados
tox
# Esperado: Passa em todas as versões
```

7️⃣ Pre-commit

```bash
# Test pre-commit hooks
pre-commit run --all-files
# Esperado: Passa sem erros
```

---

✅ Checklist Final

Após validar tudo, confirme:

- [ ] Estrutura de pastas criada (25+ arquivos)
- [ ] poetry install rodou sem erros
- [ ] python -c "from grimperium..." importa OK
- [ ] ruff check . sem critical errors
- [ ] black --check . OK
- [ ] mypy src/ OK
- [ ] pytest tests/ 100% pass
- [ ] tox passa (se você tiver múltiplas Python versions)
- [ ] pre-commit run --all-files OK

Se tudo passou ✅ → Você está pronto para o próximo batch!

---

📖 Se Algo Deu Errado

Erro: poetry install falha
Solução:
1. Delete poetry.lock
2. Delete venv: rm -rf .venv
3. poetry install novamente

Erro: ruff check com critical errors
Solução:
1. Descreva o erro para Claude Code: @debug
2. Claude Code vai tentar fixar
3. Ou rode ruff check --fix .

Erro: mypy complaining sobre types
Esperado! Tipo hints em stubs podem ter warnings. Isso é OK por enquanto.

Erro: pytest has failures
Solução:
1. Rode: pytest -v para ver qual teste falhou
2. Copie o traceback
3. Descreva para Claude Code: "Test X failed: [traceback]"

---

🎯 Próximos Passos (Após v0.1)

Depois que scaffolding está 100% validado:

1. Commit:
   ```bash
   git add .
   git commit -m "chore: v0.1 scaffolding - structure + config + stubs"
   git push origin main
   ```

2. Prepare próximo prompt para Batch 2:
   - Implementar ChemperiumLoader
   - Implementar DataFusion
   - Criar testes reais

3. Avise quando quiser o Batch 2!

---

🎬 TL;DR (Resumão Rápido)

```bash
# 1. Copie o prompt: 02_PROMPT_001_SCAFFOLDING_INICIAL.md

# 2. No Cursor: Cmd+K (ou Ctrl+K) → cole o prompt → Enter

# 3. Espere ~25 min (Claude Code vai criar tudo)

# 4. Valide:
poetry install
pytest tests/ -v
ruff check .

# 5. Se passou tudo: ✅ Sucesso!

# 6. Commit e próximo batch
git add . && git commit -m "v0.1 scaffolding"
```

---

Você está 100% pronto para começar!

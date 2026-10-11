# Dashboard local

`index.html` é a entrada do Vite. Abrir esse arquivo com `file://` não executa o
aplicativo: o navegador não transforma os módulos TypeScript/React nem resolve
o caminho `/src/main.tsx`.

Para ver dados reais, inicie a API a partir da raiz do repositório:

```bash
poetry run uvicorn grimperium.webapi.app:create_app --factory --host 127.0.0.1 --port 8000
```

Confira `http://127.0.0.1:8000/api/health` (deve retornar `status: ok`). Em
outro terminal, inicie o frontend:

```bash
cd apps/dashboard
npm run dev
```

Abra o endereço HTTP mostrado pelo Vite (normalmente
`http://localhost:5173/`). O Vite encaminha `/api` para a porta 8000. Se a API
ou um CSV estiver indisponível, a tela mostra o erro e oferece recarga.

Em um aplicativo com suporte a preview de servidores, `.claude/launch.json`
define as duas configurações: inicie `grimperium-api` e depois
`grimperium-dashboard`. As portas são fixas porque o proxy do Vite aponta para
`127.0.0.1:8000`.

## Operação

- **Progresso:** estados e qualidade CBS das linhas do `thermo_pm7.csv`,
  identificadas como protocolo legado `crest_pm7`; contagens do store
  Semi-Imperium. O botão “Atualizar dados” relê as fontes.
- **Experimental:** cobertura por fonte/família e flags, comparação CBS ×
  experimental, erro por número de átomos pesados e tabela de pares ordenada
  pelo maior erro absoluto. Busque por ID/SMILES/fonte e filtre por fonte,
  validação ou conflito. Valores `validation_only` servem apenas para
  validação, não para treino.

Todas as rotas são somente leitura. Os gráficos apresentam os pares que
correspondem aos filtros ativos.

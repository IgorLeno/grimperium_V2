# Dashboard local

`index.html` é a entrada do Vite. Abrir esse arquivo com `file://` não executa o
aplicativo: o navegador não transforma os módulos TypeScript/React nem resolve
o caminho `/src/main.tsx`.

Para desenvolvimento, inicie o frontend a partir deste diretório:

```bash
cd apps/dashboard
npm run dev
```

Abra o endereço HTTP mostrado pelo Vite (normalmente
`http://localhost:5173/`). A página atual mostra apenas o scaffold; as telas
de Progresso e Experimental são o próximo lote.

Para servir também os dados da API, inicie o backend a partir da raiz do
repositório em outro terminal:

```bash
poetry run uvicorn grimperium.webapi.app:create_app --factory --host 127.0.0.1 --port 8000
```

O Vite encaminha requisições `/api` para essa porta. O frontend funciona sem
o backend enquanto exibe apenas o scaffold.

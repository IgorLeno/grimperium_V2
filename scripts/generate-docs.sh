#!/usr/bin/env bash
set -euo pipefail

# Gera/valida a documentação do projeto via Sphinx.
#
# Uso:
#   ./scripts/generate-docs.sh
#
# Pré-requisito (uma vez por ambiente):
#   poetry install --with docs --no-interaction

show_help() {
  cat <<'EOF'
Uso: ./scripts/generate-docs.sh

Opções:
  -h, --help      Mostra esta ajuda.
EOF
}

for arg in "$@"; do
  case "$arg" in
    -h|--help) show_help; exit 0 ;;
    *) echo "Argumento desconhecido: $arg" >&2; show_help; exit 2 ;;
  esac
done

export POETRY_NO_INTERACTION=1
export PYTHONUNBUFFERED=1

if ! command -v poetry >/dev/null 2>&1; then
  echo "Erro: poetry não encontrado no PATH." >&2
  exit 127
fi

# Garante execução na raiz do repositório (útil para git hooks).
repo_root="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "$repo_root" ]]; then
  repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
cd "$repo_root"

# Validate docs dependencies are installed (Sphinx must be importable).
if ! poetry run python -c "import sphinx" >/dev/null 2>&1; then
  echo "Erro: dependências de documentação não instaladas (Sphinx não encontrado)." >&2
  echo "Execute: poetry install --with docs" >&2
  exit 1
fi

if ! poetry run sphinx-build -b html docs/source docs/build/html; then
  code=$?
  echo "Erro: sphinx-build falhou (exit code=$code)." >&2
  echo "Dica: verifique se as dependências estão instaladas com: poetry install --with docs" >&2
  exit "$code"
fi

exit 0



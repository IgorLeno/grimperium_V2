#!/usr/bin/env bash
# Instalador da célula do Grimperium (Ubuntu no WSL2, usuário comum, sem sudo).
#
# Rode de dentro da pasta extraída do pacote (pode estar em /mnt/c/...):
#
#     bash instalar.sh
#
# O que ele faz, nesta ordem:
#   1. micromamba 2.9.0 em ~/.local/bin (sha256 conferido);
#   2. ambiente grimperium-quimica (CREST, xtb, MOPAC) a partir do
#      conda-lock.yml: mesmas versões e mesmo sha256 da central;
#   3. um Python separado (ambiente grimperium-celula-python + venv em
#      ~/grimperium-celula/venv) para o programa `celula`. Fora do
#      grimperium-quimica de propósito: lá ele mudaria o environment_sha256
#      e a central recusaria o registro (409);
#   4. ~/grimperium-celula/ com config.json (chmod 600) e iniciar.sh;
#   5. confere os programas contra os da central e se a central responde;
#   6. no WSL, deixa o script do Windows (.wslconfig, Agendador, energia) na
#      pasta do seu usuário do Windows.
#
# Rodar de novo é seguro: refaz os ambientes e troca o token pelo do pacote
# novo (pare a célula antes).
set -euo pipefail

PACOTE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CELULA="${GRIMPERIUM_CELULA:-$HOME/grimperium-celula}"
MAMBA_ROOT="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
BIN="$HOME/.local/bin"
MICROMAMBA="$BIN/micromamba"
MICROMAMBA_URL="https://github.com/mamba-org/micromamba-releases/releases/download/2.9.0-0/micromamba-linux-64"
MICROMAMBA_SHA256="366cd9cd8be14df1ab8ed50352a82111082a36686b2d389fdb79a92c3fafb3e3"
AMBIENTE="grimperium-quimica"
AMBIENTE_PYTHON="grimperium-celula-python"

passo() { printf '\n==> %s\n' "$*"; }
falha() { printf '\nERRO: %s\n' "$*" >&2; exit 1; }
aviso() { printf 'AVISO: %s\n' "$*" >&2; }

baixar() {  # baixar URL DESTINO
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL --retry 3 -o "$2" "$1"
    elif command -v wget >/dev/null 2>&1; then
        wget -q -O "$2" "$1"
    else
        falha "nem curl nem wget instalados; rode: sudo apt install curl"
    fi
}

mm() { "$MICROMAMBA" -r "$MAMBA_ROOT" "$@"; }

# --- verificações -----------------------------------------------------------

[ "$(id -u)" -ne 0 ] || falha "não rode como root nem com sudo; use seu usuário"
[ "$(uname -m)" = "x86_64" ] || falha "só computadores x86_64 (Intel/AMD)"
for arquivo in pacote.env pacote.json config.json conda-lock.yml \
    requirements-celula.txt constraints-central.txt; do
    [ -f "$PACOTE/$arquivo" ] || falha "pacote incompleto: falta $arquivo"
done
# shellcheck source=/dev/null
. "$PACOTE/pacote.env"
[ -f "$PACOTE/$WHEEL" ] || falha "pacote incompleto: falta $WHEEL"
case "$CELULA/" in
    /mnt/*) falha "a pasta da célula ($CELULA) não pode ficar num disco do Windows" ;;
esac
case "$MAMBA_ROOT/" in
    /mnt/*) falha "os ambientes ($MAMBA_ROOT) não podem ficar num disco do Windows" ;;
esac
command -v flock >/dev/null 2>&1 || falha "falta o comando flock (pacote util-linux)"
command -v setsid >/dev/null 2>&1 || falha "falta o comando setsid (pacote util-linux)"
mkdir -p "$CELULA" "$BIN"
chmod 700 "$CELULA"
if ! flock -n "$CELULA/celula.lock" true; then
    falha "a célula está rodando; pare-a antes (Ctrl+C na janela dela, ou reinicie o WSL)"
fi

echo "Célula do Grimperium: $CELL_NAME"
echo "Central: $CENTRAL_URL"
echo "Pasta:   $CELULA"

# --- 1. micromamba ----------------------------------------------------------

passo "1/6 micromamba"
if [ -x "$MICROMAMBA" ] && echo "$MICROMAMBA_SHA256  $MICROMAMBA" | sha256sum -c --status; then
    echo "já instalado em $MICROMAMBA"
else
    temporario="$(mktemp)"
    baixar "$MICROMAMBA_URL" "$temporario"
    if ! echo "$MICROMAMBA_SHA256  $temporario" | sha256sum -c --status; then
        rm -f "$temporario"
        falha "o micromamba baixado não tem o sha256 esperado; tente de novo mais tarde"
    fi
    install -m 755 "$temporario" "$MICROMAMBA"
    rm -f "$temporario"
    echo "instalado em $MICROMAMBA"
fi

# --- 2. ambiente dos programas ----------------------------------------------

passo "2/6 ambiente $AMBIENTE (CREST, xtb, MOPAC; pode levar alguns minutos)"
if [ -d "$MAMBA_ROOT/envs/$AMBIENTE" ]; then
    mm env remove -y -n "$AMBIENTE" >/dev/null
fi
mm create -y -n "$AMBIENTE" -f "$PACOTE/conda-lock.yml"

# --- 3. Python da célula ----------------------------------------------------

passo "3/6 Python $PYTHON_VERSION da célula (separado dos programas)"
if [ -d "$MAMBA_ROOT/envs/$AMBIENTE_PYTHON" ]; then
    mm env remove -y -n "$AMBIENTE_PYTHON" >/dev/null
fi
mm create -y -n "$AMBIENTE_PYTHON" -c conda-forge "python=$PYTHON_VERSION"
rm -rf "$CELULA/venv"
"$MAMBA_ROOT/envs/$AMBIENTE_PYTHON/bin/python" -m venv "$CELULA/venv"
PY="$CELULA/venv/bin/python"
"$PY" -m pip install --disable-pip-version-check --no-input \
    -r "$PACOTE/requirements-celula.txt" -c "$PACOTE/constraints-central.txt"
"$PY" -m pip install --disable-pip-version-check --no-input \
    --no-deps --force-reinstall "$PACOTE/$WHEEL"

# --- 4. pasta da célula -----------------------------------------------------

passo "4/6 pasta $CELULA"
# Mantém ajustes locais (workers, crest_threads, ...); central e token vêm do
# pacote. O arquivo nasce 600: o token nunca fica legível por outros.
umask 077
"$PY" - "$PACOTE/config.json" "$CELULA/config.json" <<'PYTHON'
import json, os, sys
novo = json.load(open(sys.argv[1], encoding="utf-8"))
destino = sys.argv[2]
atual = {}
if os.path.exists(destino):
    try:
        atual = json.load(open(destino, encoding="utf-8"))
    except ValueError:
        atual = {}
atual.update(novo)
temporario = destino + ".tmp"
with open(os.open(temporario, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600), "w", encoding="utf-8") as arquivo:
    json.dump(atual, arquivo, indent=2)
    arquivo.write("\n")
os.chmod(temporario, 0o600)
os.replace(temporario, destino)
PYTHON
umask 022
install -m 644 "$PACOTE/pacote.json" "$CELULA/pacote.json"
cat > "$CELULA/iniciar.sh" <<INICIAR
#!/usr/bin/env bash
# Chamado pela tarefa do Agendador do Windows. setsid -w: a célula lidera o
# próprio grupo de processos (o segundo sinal de parada leva o CREST junto) e
# o wsl.exe fica aberto enquanto ela roda, o que mantém o WSL ligado.
exec setsid -w "$CELULA/venv/bin/celula" --config "$CELULA/config.json" start
INICIAR
chmod 755 "$CELULA/iniciar.sh"
ln -sf "$CELULA/venv/bin/celula" "$BIN/celula"
echo "config.json (600), pacote.json e iniciar.sh prontos; comando: $BIN/celula"

# --- 5. conferências --------------------------------------------------------

passo "5/6 programas iguais aos da central?"
if ! "$PY" - "$CELULA/pacote.json" <<'PYTHON'
import json, sys
from semi_imperium.training.__main__ import _executables
from semi_imperium.training.program_info import probe_programs, program_mismatches
central = json.load(open(sys.argv[1], encoding="utf-8"))["central_programs"]
celula = probe_programs(_executables(None))
diferencas = program_mismatches(central, celula)
for chave in ("crest", "xtb", "mopac"):
    print(f"  {chave} {celula.get(chave)}")
print(f"  environment_sha256 {celula.get('environment_sha256', '')[:16]}…")
if diferencas:
    print("diferente da central: " + "; ".join(diferencas), file=sys.stderr)
    sys.exit(1)
print("iguais aos da central")
PYTHON
then
    falha "os programas não batem com os da central: o registro seria recusado (409)"
fi

passo "6/6 a central responde?"
if command -v curl >/dev/null 2>&1; then
    codigo="$(curl -s -o /dev/null -w '%{http_code}' --max-time 10 "$CENTRAL_URL/api/cell/register" || true)"
    if [ "$codigo" = "000" ] || [ -z "$codigo" ]; then
        aviso "sem resposta de $CENTRAL_URL. Confira se o Tailscale está ligado nos dois computadores e se a API das células está rodando na central."
    else
        echo "central respondeu (HTTP $codigo)"
    fi
else
    aviso "sem curl: não deu para testar a conexão com a central"
fi

# --- Windows ----------------------------------------------------------------

SCRIPT_WINDOWS=""
if [ -n "${WSL_DISTRO_NAME:-}" ] && command -v cmd.exe >/dev/null 2>&1 \
    && command -v wslpath >/dev/null 2>&1; then
    perfil="$(cd /mnt/c 2>/dev/null && cmd.exe /c 'echo %USERPROFILE%' 2>/dev/null | tr -d '\r' || true)"
    if [ -n "$perfil" ]; then
        SCRIPT_WINDOWS="$(wslpath "$perfil")/grimperium-celula-windows.ps1"
        # BOM: o PowerShell 5 do Windows lê script sem BOM como ANSI e
        # estraga os acentos.
        {
            printf '\xef\xbb\xbf'
            sed -e "s|@DISTRO@|$WSL_DISTRO_NAME|g" \
                -e "s|@LINUX_USER@|$(id -un)|g" \
                -e "s|@INICIAR@|$CELULA/iniciar.sh|g" \
                "$PACOTE/windows/configurar-windows.ps1"
        } > "$SCRIPT_WINDOWS"
    fi
fi

printf '\nInstalação concluída.\n\n'
echo "Próximos passos:"
if [ -n "$SCRIPT_WINDOWS" ]; then
    echo "  1. No Windows, abra o PowerShell (sem ser administrador) e rode:"
    echo "       powershell -ExecutionPolicy Bypass -File \"$perfil\\grimperium-celula-windows.ps1\""
    echo "     Ele ajusta o .wslconfig, cria a tarefa no Agendador e tira a suspensão."
    echo "     Depois reinicie o Windows: a célula sobe sozinha no logon."
else
    echo "  1. Isto não parece o WSL: a parte do Windows foi pulada."
    echo "     Para iniciar à mão: $CELULA/iniciar.sh"
fi
echo "  2. A célula faz sozinha o teste de admissão (5 referências) e então recebe a fila."
echo "     Acompanhe na central: dashboard, Calcular PM7 > Células."
echo "  3. Acompanhe aqui com: celula status   (log: $CELULA/celula.log)"
echo
echo "O pacote extraído ($PACOTE) e o zip contêm o token da célula: apague os dois."

#!/usr/bin/env bash
# Gera dist/OpenDTA.app e dist/OpenDTA-<versão>.dmg no macOS.
#
#   tools/build_macos.sh            app + dmg
#   tools/build_macos.sh --no-dmg   só o app
#
# Usa um ambiente virtual próprio (.venv-build), para o pacote não levar
# bibliotecas de desenvolvimento. Assinatura ad hoc: o app roda na máquina
# em que foi gerado; para distribuir a outras pessoas seria preciso uma
# conta de desenvolvedor da Apple (assinatura e notarização).
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$(pwd)"

if [[ "$(uname)" != "Darwin" ]]; then
    echo "Este script é para o macOS. Em outros sistemas use: pyinstaller packaging/opendta.spec" >&2
    exit 1
fi

PY="${PYTHON:-python3}"
"$PY" - <<'PYCHECK'
import sys
if sys.version_info < (3, 11):
    sys.exit("O OpenDTA precisa do Python 3.11 ou mais recente.")
PYCHECK

if [[ ! -d .venv-build ]]; then
    "$PY" -m venv .venv-build
fi
source .venv-build/bin/activate
python -m pip install --upgrade pip >/dev/null
python -m pip install . "pyinstaller>=6.3" pyreadstat

VERSION="$(python -c 'import opendta; print(opendta.__version__)')"
rm -rf build dist/OpenDTA dist/OpenDTA.app
pyinstaller packaging/opendta.spec --noconfirm --clean

# assinatura ad hoc (o macOS exige alguma assinatura em Apple Silicon)
codesign --force --deep --sign - dist/OpenDTA.app

echo
echo "Aplicativo: $ROOT/dist/OpenDTA.app"

if [[ "${1:-}" != "--no-dmg" ]]; then
    # DMG_NAME permite nomes fixos no GitHub (ex.: OpenDTA-macos-apple-silicon)
    DMG="dist/${DMG_NAME:-OpenDTA-$VERSION}.dmg"
    STAGE="$(mktemp -d)"
    cp -R dist/OpenDTA.app "$STAGE/"
    ln -s /Applications "$STAGE/Applications"
    cp LICENSE "$STAGE/LICENSE.txt"
    rm -f "$DMG"
    # o hdiutil às vezes falha com "Resource busy" nas máquinas do GitHub
    for tentativa in 1 2 3 4 5; do
        if hdiutil create -volname "OpenDTA $VERSION" -srcfolder "$STAGE" -ov -format UDZO "$DMG" >/dev/null; then
            break
        fi
        echo "hdiutil falhou (tentativa $tentativa); tentando de novo..." >&2
        sleep 5
    done
    rm -rf "$STAGE"
    [[ -f "$DMG" ]] || { echo "Não foi possível criar o .dmg" >&2; exit 1; }
    echo "Imagem de disco: $ROOT/$DMG"
fi

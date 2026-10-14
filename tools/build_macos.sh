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
    DMG="dist/OpenDTA-$VERSION.dmg"
    STAGE="$(mktemp -d)"
    cp -R dist/OpenDTA.app "$STAGE/"
    ln -s /Applications "$STAGE/Applications"
    rm -f "$DMG"
    hdiutil create -volname "OpenDTA $VERSION" -srcfolder "$STAGE" -ov -format UDZO "$DMG" >/dev/null
    rm -rf "$STAGE"
    echo "Imagem de disco: $ROOT/$DMG"
fi

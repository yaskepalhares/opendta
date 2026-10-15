#!/usr/bin/env bash
# Builds dist/OpenDTA.app and dist/OpenDTA-<version>.dmg on macOS.
#
#   tools/build_macos.sh            app + dmg
#   tools/build_macos.sh --no-dmg   the app only
#
# Uses its own virtual environment (.venv-build), so that the bundle does not
# carry development libraries. Ad hoc signature: the app runs on the machine
# where it was built; distributing it to other people would need an Apple
# developer account (signing and notarization).
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$(pwd)"

if [[ "$(uname)" != "Darwin" ]]; then
    echo "This script is for macOS. On other systems use: pyinstaller packaging/opendta.spec" >&2
    exit 1
fi

PY="${PYTHON:-python3}"
"$PY" - <<'PYCHECK'
import sys
if sys.version_info < (3, 11):
    sys.exit("OpenDTA needs Python 3.11 or later.")
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

# ad hoc signature (macOS requires some signature on Apple Silicon)
codesign --force --deep --sign - dist/OpenDTA.app

echo
echo "Application: $ROOT/dist/OpenDTA.app"

if [[ "${1:-}" != "--no-dmg" ]]; then
    # DMG_NAME allows fixed names on GitHub (e.g. OpenDTA-macos-apple-silicon)
    DMG="dist/${DMG_NAME:-OpenDTA-$VERSION}.dmg"
    STAGE="$(mktemp -d)"
    cp -R dist/OpenDTA.app "$STAGE/"
    ln -s /Applications "$STAGE/Applications"
    cp LICENSE "$STAGE/LICENSE.txt"
    rm -f "$DMG"
    # hdiutil sometimes fails with "Resource busy" on GitHub machines
    for attempt in 1 2 3 4 5; do
        if hdiutil create -volname "OpenDTA $VERSION" -srcfolder "$STAGE" -ov -format UDZO "$DMG" >/dev/null; then
            break
        fi
        echo "hdiutil failed (attempt $attempt); trying again..." >&2
        sleep 5
    done
    rm -rf "$STAGE"
    [[ -f "$DMG" ]] || { echo "Could not create the .dmg" >&2; exit 1; }
    echo "Disk image: $ROOT/$DMG"
fi

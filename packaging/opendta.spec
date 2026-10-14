# Especificação do PyInstaller para o OpenDTA.
#
#   pyinstaller packaging/opendta.spec --noconfirm
#
# No macOS gera dist/OpenDTA.app (pasta única, assinatura ad hoc); no Linux e
# no Windows, a pasta dist/OpenDTA com o executável. tools/build_macos.sh faz
# o processo inteiro e ainda monta um .dmg.
# -*- mode: python -*-

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

ROOT = Path(SPECPATH).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from opendta import __version__  # noqa: E402

datas = collect_data_files("opendta", includes=["help/*.sthlp", "gui/icons/*.svg", "gui/appicon/*.png"])
hidden = collect_submodules("opendta")

# licenças: a do OpenDTA, o aviso de terceiros e os metadados (com o texto das
# licenças) dos pacotes levados no app
datas += [(str(ROOT / "LICENSE"), "."), (str(ROOT / "packaging" / "THIRD_PARTY_NOTICES.md"), ".")]
for pkg in ("PySide6", "shiboken6", "numpy", "scipy", "matplotlib", "openpyxl", "pillow",
            "pyreadstat", "pandas", "fonttools", "kiwisolver", "contourpy", "cycler", "pyparsing",
            "python-dateutil", "packaging", "et-xmlfile"):
    try:
        datas += copy_metadata(pkg)
    except Exception:  # noqa: BLE001 - pacote opcional ausente no ambiente de build
        pass

a = Analysis(
    [str(ROOT / "packaging" / "launcher.py")],
    pathex=[str(ROOT / "src")],
    datas=datas,
    hiddenimports=hidden + ["matplotlib.backends.backend_qtagg", "matplotlib.backends.backend_agg",
                            "matplotlib.backends.backend_pdf", "matplotlib.backends.backend_svg",
                            "matplotlib.backends.backend_ps", "scipy.special._cdflib"],
    excludes=["tkinter", "PyQt5", "PyQt6", "IPython", "pandas.tests", "pytest"],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="OpenDTA",
    console=False,               # sem janela de terminal
    argv_emulation=False,        # os arquivos chegam pelo QFileOpenEvent (gui/app.py)
    icon=str(ROOT / "assets/icons/opendta/windows/opendta.ico")
    if sys.platform == "win32" else None,
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="OpenDTA")

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="OpenDTA.app",
        icon=str(ROOT / "assets/icons/opendta/macos/opendta.icns"),
        bundle_identifier="io.github.yaskepalhares.opendta",
        version=__version__,
        info_plist={
            "CFBundleName": "OpenDTA",
            "CFBundleDisplayName": "OpenDTA",
            "CFBundleShortVersionString": __version__,
            "NSHighResolutionCapable": True,
            "NSRequiresAquaSystemAppearance": False,
            "LSMinimumSystemVersion": "12.0",      # mínimo do Qt 6.8
            "NSHumanReadableCopyright": "© 2026 Yaske Palhares. MIT License. "
                                        "Not affiliated with StataCorp.",
            # "Abrir com" para do-files e bases .dta; o OpenDTA não se torna o
            # aplicativo padrão desses arquivos (LSHandlerRank = Alternate)
            "CFBundleDocumentTypes": [
                {"CFBundleTypeName": "Do-file", "CFBundleTypeRole": "Editor",
                 "LSHandlerRank": "Alternate", "CFBundleTypeExtensions": ["do", "ado"]},
                {"CFBundleTypeName": "Dataset (.dta)", "CFBundleTypeRole": "Editor",
                 "LSHandlerRank": "Alternate", "CFBundleTypeExtensions": ["dta"]},
            ],
        },
    )

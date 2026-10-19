# PyInstaller spec for OpenDTA.
#
#   pyinstaller packaging/opendta.spec --noconfirm
#
# On macOS it builds dist/OpenDTA.app (one folder, ad hoc signature); on Linux and
# Windows, the dist/OpenDTA folder with the executable. tools/build_macos.sh runs
# the whole process and also builds a .dmg.
# -*- mode: python -*-

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

ROOT = Path(SPECPATH).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from opendta import __version__  # noqa: E402

datas = collect_data_files("opendta", includes=["help/*.sthlp", "gui/icons/*.svg", "gui/appicon/*.png", "data/*.dta", "data/icd/*"])
# certifi's CA bundle (cacert.pem), used by the update check
datas += collect_data_files("certifi")
hidden = collect_submodules("opendta")

# licenses: OpenDTA's own, the third-party notice and the metadata (with the license texts) of the
# packages bundled in the app
datas += [(str(ROOT / "LICENSE"), "."), (str(ROOT / "packaging" / "THIRD_PARTY_NOTICES.md"), ".")]
for pkg in ("PySide6", "shiboken6", "numpy", "scipy", "matplotlib", "openpyxl", "pillow",
            "pyreadstat", "pandas", "fonttools", "kiwisolver", "contourpy", "cycler", "pyparsing",
            "python-dateutil", "packaging", "et-xmlfile", "certifi"):
    try:
        datas += copy_metadata(pkg)
    except Exception:  # noqa: BLE001 - optional package missing from the build environment
        pass

a = Analysis(
    [str(ROOT / "packaging" / "launcher.py")],
    pathex=[str(ROOT / "src")],
    datas=datas,
    hiddenimports=hidden + ["matplotlib.backends.backend_qtagg", "matplotlib.backends.backend_agg",
                            "matplotlib.backends.backend_pdf", "matplotlib.backends.backend_svg",
                            "matplotlib.backends.backend_ps", "scipy.special._cdflib",
                            # imported through opendta.lazy (names in strings, invisible to the
                            # analysis)
                            "scipy.stats", "scipy.special", "scipy.linalg", "matplotlib.figure",
                            "matplotlib.ticker"],
    excludes=["tkinter", "PyQt5", "PyQt6", "IPython", "pandas.tests", "pytest"],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="OpenDTA",
    console=False,               # no terminal window
    argv_emulation=False,        # files arrive through QFileOpenEvent (gui/app.py)
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
            "LSMinimumSystemVersion": "12.0",      # minimum for Qt 6.8
            "NSHumanReadableCopyright": "© 2026 OpenDTA contributors. MIT License. "
                                        "Not affiliated with StataCorp.",
            # "Open With" for do-files and .dta datasets; OpenDTA does not become the default
            # application for these files (LSHandlerRank = Alternate)
            "CFBundleDocumentTypes": [
                {"CFBundleTypeName": "Do-file", "CFBundleTypeRole": "Editor",
                 "LSHandlerRank": "Alternate", "CFBundleTypeExtensions": ["do", "ado"]},
                {"CFBundleTypeName": "Dataset (.dta)", "CFBundleTypeRole": "Editor",
                 "LSHandlerRank": "Alternate", "CFBundleTypeExtensions": ["dta"]},
            ],
        },
    )

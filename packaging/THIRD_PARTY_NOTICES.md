# Third-party software

The OpenDTA application bundle (`OpenDTA.app` and the Windows build) contains
the following open-source components. OpenDTA itself is released under the MIT
License (see `LICENSE`). The full license text of each Python package is
included in the bundle, in its `*.dist-info` folder inside `Contents/Resources`
(macOS) or `_internal` (Windows/Linux).

| Component | License | Source |
|---|---|---|
| Python | PSF License | https://www.python.org |
| Qt 6 and PySide6 / Shiboken6 | LGPL-3.0 | https://www.qt.io, https://wiki.qt.io/Qt_for_Python |
| NumPy | BSD-3-Clause | https://numpy.org |
| SciPy | BSD-3-Clause | https://scipy.org |
| Matplotlib | Matplotlib License (PSF-based) | https://matplotlib.org |
| Pillow | MIT-CMU (HPND) | https://python-pillow.org |
| contourpy, cycler, kiwisolver | BSD-3-Clause | https://pypi.org |
| fonttools | MIT | https://github.com/fonttools/fonttools |
| pyparsing | MIT | https://github.com/pyparsing/pyparsing |
| python-dateutil | Apache-2.0 / BSD-3-Clause | https://github.com/dateutil/dateutil |
| packaging | Apache-2.0 / BSD-2-Clause | https://github.com/pypa/packaging |
| openpyxl, et-xmlfile | MIT | https://openpyxl.readthedocs.io |
| pyreadstat | Apache-2.0 | https://github.com/Roche/pyreadstat |
| ReadStat | MIT | https://github.com/WizardMac/ReadStat |
| pandas | BSD-3-Clause | https://pandas.pydata.org |
| PyInstaller bootloader | GPL-2.0 with a bootloader exception | https://pyinstaller.org |

## Qt and the LGPL

Qt and PySide6 are used under the GNU Lesser General Public License v3. They are
shipped as separate, unmodified shared libraries inside the application folder
(`Contents/Frameworks` on macOS), so they can be replaced by other compatible
builds. The source code of Qt is available at https://download.qt.io and the
source of PySide6 at https://code.qt.io/cgit/pyside/pyside-setup.git. The LGPL
text is at https://www.gnu.org/licenses/lgpl-3.0.html.

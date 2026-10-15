# Builds the OpenDTA Windows installer (64-bit, Intel/AMD) and a portable .zip.
#
#   powershell -ExecutionPolicy Bypass -File tools\build_windows.ps1
#
# Produces dist\OpenDTA-windows-x64-setup.exe (Inno Setup) and
# dist\OpenDTA-windows-x64-portable.zip. Uses its own virtual environment
# (.venv-build). The installer is not signed: Windows SmartScreen asks for
# confirmation on first run ("More info" -> "Run anyway").

$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..")
$Root = Get-Location

$Python = if ($env:PYTHON) { $env:PYTHON } else { "python" }
& $Python -c "import sys, platform; sys.exit(0 if sys.version_info >= (3, 11) and platform.machine().lower() in ('amd64', 'x86_64') else 1)"
if ($LASTEXITCODE -ne 0) { throw "OpenDTA needs 64-bit (x64) Python 3.11 or later." }

if (-not (Test-Path ".venv-build")) { & $Python -m venv .venv-build }
$Py = Join-Path $Root ".venv-build\Scripts\python.exe"
& $Py -m pip install --upgrade pip | Out-Null
& $Py -m pip install . "pyinstaller>=6.3" pyreadstat
if ($LASTEXITCODE -ne 0) { throw "pip install failed" }

$Version = & $Py -c "import opendta; print(opendta.__version__)"
Remove-Item -Recurse -Force build, dist\OpenDTA -ErrorAction SilentlyContinue
& $Py -m PyInstaller packaging\opendta.spec --noconfirm --clean
if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed" }

# portable version: the PyInstaller folder, zipped
$Zip = "dist\OpenDTA-windows-x64-portable.zip"
Remove-Item $Zip -ErrorAction SilentlyContinue
Copy-Item LICENSE dist\OpenDTA\LICENSE.txt
Compress-Archive -Path dist\OpenDTA -DestinationPath $Zip
Write-Host "Portable: $Root\$Zip"

# installer: Inno Setup 6 (installed with Chocolatey if missing)
$Iscc = Get-Command iscc.exe -ErrorAction SilentlyContinue
if (-not $Iscc) {
    $Candidates = @("${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe", "$env:ProgramFiles\Inno Setup 6\ISCC.exe",
                    "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe")
    $Found = $Candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
    if (-not $Found) {
        if (Get-Command choco -ErrorAction SilentlyContinue) {
            choco install innosetup -y --no-progress | Out-Null
            $Found = $Candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
        }
    }
    if (-not $Found) { throw "Inno Setup 6 not found: install it from https://jrsoftware.org/isdl.php" }
    $Iscc = $Found
} else { $Iscc = $Iscc.Source }
& $Iscc "/DMyAppVersion=$Version" packaging\opendta.iss
if ($LASTEXITCODE -ne 0) { throw "Inno Setup failed" }
Write-Host "Installer: $Root\dist\OpenDTA-windows-x64-setup.exe"

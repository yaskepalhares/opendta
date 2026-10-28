"""Update check and installation from the GitHub Releases of OpenDTA.

The latest release is read from the GitHub REST API; its version is
compared with opendta.__version__ and, when newer, the asset for this
installation is downloaded (checked against its size and, when GitHub
gives it, its SHA-256 digest) and installed after OpenDTA closes:

- Windows, installed: the Inno Setup installer runs silently over the
  current installation (same AppId) and opens OpenDTA again (/RELAUNCH=1).
- Windows, portable: the new .zip is unpacked and copied over the folder
  by a small batch script that waits for OpenDTA to exit.
- macOS: the .app is copied out of the .dmg and a shell script replaces the
  running bundle after OpenDTA exits, then opens it. When the folder is not
  writable, the .dmg is opened in the Finder instead.
- Running from source: only the notice (update with git and pip).
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import ssl
import subprocess
import sys
import tempfile
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from . import __version__

REPO = "yaskepalhares/opendta"
API_LATEST = f"https://api.github.com/repos/{REPO}/releases/latest"
RELEASES_PAGE = f"https://github.com/{REPO}/releases/latest"

ASSETS = {
    "macos-arm64": "OpenDTA-macos-apple-silicon.dmg",
    "macos-x86_64": "OpenDTA-macos-intel.dmg",
    "windows-installer": "OpenDTA-windows-x64-setup.exe",
    "windows-portable": "OpenDTA-windows-x64-portable.zip",
}


class UpdateError(Exception):
    pass


@dataclass
class Asset:
    name: str
    url: str
    size: int = 0
    digest: str = ""          # "sha256:<hex>" when GitHub provides it


@dataclass
class Release:
    tag: str
    version: str
    notes: str = ""
    page: str = RELEASES_PAGE
    assets: dict[str, Asset] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# versions
# ---------------------------------------------------------------------------

def parse_version(text: str) -> tuple:
    """'v0.2.1' -> (0, 2, 1, 1); a pre-release ('0.3.0-beta2', '0.3.0rc1') sorts
    before the release: (0, 3, 0, 0, 'beta', 2)."""
    t = text.strip().lstrip("vV")
    m = re.match(r"(\d+(?:\.\d+)*)(.*)$", t)
    if not m:
        return (0,)
    nums = [int(x) for x in m.group(1).split(".")]
    while len(nums) < 3:
        nums.append(0)
    rest = m.group(2).strip(".-_+ ")
    if not rest:
        return tuple(nums) + (1,)
    pm = re.match(r"([A-Za-z]+)\.?(\d*)", rest)
    label, n = (pm.group(1).lower(), int(pm.group(2) or 0)) if pm else (rest, 0)
    return tuple(nums) + (0, label, n)


def is_newer(remote: str, local: str = __version__) -> bool:
    return parse_version(remote) > parse_version(local)


# ---------------------------------------------------------------------------
# where and how OpenDTA is installed
# ---------------------------------------------------------------------------

def install_kind() -> str:
    """'macos-arm64', 'macos-x86_64', 'windows-installer', 'windows-portable' or 'source'."""
    if not getattr(sys, "frozen", False):
        return "source"
    if sys.platform == "darwin":
        return "macos-arm64" if platform.machine() == "arm64" else "macos-x86_64"
    if sys.platform == "win32":
        exe_dir = Path(sys.executable).resolve().parent
        return "windows-installer" if any(exe_dir.glob("unins*.exe")) else "windows-portable"
    return "source"


def app_location() -> Path | None:
    """The .app bundle (macOS) or the program folder (Windows) being run."""
    exe = Path(sys.executable).resolve()
    if sys.platform == "darwin":
        for p in exe.parents:
            if p.suffix == ".app":
                return p
        return None
    return exe.parent


# ---------------------------------------------------------------------------
# GitHub
# ---------------------------------------------------------------------------

def _request(url: str) -> urllib.request.Request:
    return urllib.request.Request(url, headers={
        "Accept": "application/vnd.github+json",
        "User-Agent": f"OpenDTA/{__version__}",
        "X-GitHub-Api-Version": "2022-11-28",
    })


# The release notes hold both languages: English first, then a "## Português" section. Each
# part starts with (or is preceded by) a line that links to the other language.
_PT_HEADING = re.compile(r"^#{1,6}[ \t]*Portugu[eê]s[ \t]*$", re.MULTILINE)
_LANG_LINKS = re.compile(r"^.*English.*Portugu[eê]s.*$\n?", re.MULTILINE)
_RULE_AT_END = re.compile(r"\n-{3,}[ \t]*\Z")


def notes_for(body: str, lang: str) -> str:
    """The part of the release notes in `lang` ("pt" or anything else for English). Notes without
    a Portuguese section come back whole, minus the language links."""
    m = _PT_HEADING.search(body)
    if m is None:
        text = body
    elif lang == "pt":
        text = body[m.end():]
    else:
        text = _RULE_AT_END.sub("", body[:m.start()].rstrip())
    return _LANG_LINKS.sub("", text).strip()


def parse_release(data: dict) -> Release:
    tag = str(data.get("tag_name", ""))
    rel = Release(tag=tag, version=tag.lstrip("vV"), notes=str(data.get("body") or ""),
                  page=str(data.get("html_url") or RELEASES_PAGE))
    for a in data.get("assets") or []:
        rel.assets[a.get("name", "")] = Asset(a.get("name", ""), a.get("browser_download_url", ""),
                                              int(a.get("size") or 0), str(a.get("digest") or ""))
    return rel


def ssl_context() -> ssl.SSLContext:
    """Certificate check with certifi's CA bundle. The Python inside the packaged application
    looks for the certificates where the build machine kept them, which do not exist on the
    user's computer ("certificate verify failed"); the system store is the fallback when certifi
    is missing."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except (ImportError, OSError):
        return ssl.create_default_context()


def _urlopen(req, timeout: float):
    return urllib.request.urlopen(req, timeout=timeout, context=ssl_context())


def fetch_latest(timeout: float = 10.0, opener: Callable | None = None) -> Release:
    """The latest published release (drafts and pre-releases are not returned by GitHub)."""
    open_ = opener or _urlopen
    try:
        with open_(_request(API_LATEST), timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as e:      # noqa: BLE001 - network, HTTP or JSON errors
        raise UpdateError(str(e)) from e
    rel = parse_release(data)
    if not rel.version:
        raise UpdateError("the release has no version tag")
    return rel


def asset_for(rel: Release, kind: str | None = None) -> Asset | None:
    kind = kind or install_kind()
    name = ASSETS.get(kind)
    return rel.assets.get(name) if name else None


# ---------------------------------------------------------------------------
# download
# ---------------------------------------------------------------------------

def download(asset: Asset, dest_dir: Path, progress: Callable[[int, int], None] | None = None,
             cancelled: Callable[[], bool] | None = None, opener: Callable | None = None) -> Path:
    """Downloads the asset to dest_dir and checks its size and SHA-256 digest."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / asset.name
    part = target.with_suffix(target.suffix + ".part")
    open_ = opener or _urlopen
    h = hashlib.sha256()
    done = 0
    try:
        req = urllib.request.Request(asset.url, headers={"User-Agent": f"OpenDTA/{__version__}"}) \
            if asset.url.startswith("http") else asset.url
        with open_(req, timeout=60) as resp, open(part, "wb") as f:
            headers = getattr(resp, "headers", None)
            total = asset.size or int((headers.get("Content-Length") if headers is not None else 0) or 0)
            while True:
                if cancelled is not None and cancelled():
                    raise UpdateError("cancelled")
                chunk = resp.read(1 << 16)
                if not chunk:
                    break
                f.write(chunk)
                h.update(chunk)
                done += len(chunk)
                if progress is not None:
                    progress(done, total)
    except UpdateError:
        part.unlink(missing_ok=True)
        raise
    except Exception as e:      # noqa: BLE001
        part.unlink(missing_ok=True)
        raise UpdateError(str(e)) from e
    if asset.size and done != asset.size:
        part.unlink(missing_ok=True)
        raise UpdateError(f"incomplete download ({done} of {asset.size} bytes)")
    if asset.digest.startswith("sha256:") and h.hexdigest() != asset.digest.split(":", 1)[1].lower():
        part.unlink(missing_ok=True)
        raise UpdateError("the downloaded file does not match its SHA-256 digest")
    os.replace(part, target)
    return target


# ---------------------------------------------------------------------------
# installation (runs after OpenDTA exits)
# ---------------------------------------------------------------------------

def plan_install(kind: str, file: Path) -> Callable[[], None]:
    """Prepares the installation and returns the function that starts it; call it right
    before OpenDTA exits. Raises UpdateError when this installation cannot be updated
    in place (the caller then opens the downloaded file or the releases page)."""
    pid = os.getpid()
    if kind == "windows-installer":
        def start():
            flags = getattr(subprocess, "DETACHED_PROCESS", 0) | getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
            subprocess.Popen([str(file), "/SILENT", "/SUPPRESSMSGBOXES", "/NORESTART", "/RELAUNCH=1"],
                             creationflags=flags, close_fds=True)
        return start
    if kind == "windows-portable":
        target = app_location()
        if target is None or not os.access(target, os.W_OK):
            raise UpdateError("the program folder is not writable")
        stage = Path(tempfile.mkdtemp(prefix="opendta-update-"))
        shutil.unpack_archive(str(file), str(stage))
        inner = next((p for p in stage.iterdir() if p.is_dir() and (p / "OpenDTA.exe").exists()), stage)
        script = stage / "update.cmd"
        script.write_text(
            "@echo off\r\n"
            f":wait\r\ntasklist /FI \"PID eq {pid}\" | find \"{pid}\" >nul && (timeout /t 1 /nobreak >nul & goto wait)\r\n"
            f"robocopy \"{inner}\" \"{target}\" /E /NFL /NDL /NJH /NJS /NP >nul\r\n"
            f"start \"\" \"{target / 'OpenDTA.exe'}\"\r\n", encoding="utf-8")

        def start():
            flags = getattr(subprocess, "DETACHED_PROCESS", 0) | getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
            subprocess.Popen(["cmd", "/c", str(script)], creationflags=flags, close_fds=True)
        return start
    if kind.startswith("macos"):
        target = app_location()
        if target is None or not os.access(target.parent, os.W_OK) or not os.access(target, os.W_OK):
            raise UpdateError("the Applications folder is not writable")
        stage = Path(tempfile.mkdtemp(prefix="opendta-update-"))
        mount = stage / "mnt"
        mount.mkdir()
        r = subprocess.run(["hdiutil", "attach", "-nobrowse", "-readonly", "-noautoopen", "-mountpoint",
                            str(mount), str(file)], capture_output=True, text=True)
        if r.returncode != 0:
            raise UpdateError(f"could not open the disk image: {r.stderr.strip()}")
        try:
            src = next(mount.glob("*.app"), None)
            if src is None:
                raise UpdateError("the disk image has no application")
            new_app = stage / src.name
            r = subprocess.run(["ditto", str(src), str(new_app)], capture_output=True, text=True)
            if r.returncode != 0:
                raise UpdateError(f"could not copy the application: {r.stderr.strip()}")
        finally:
            subprocess.run(["hdiutil", "detach", "-quiet", str(mount)], capture_output=True)
        old = target.with_name(target.name + ".old")
        script = stage / "update.sh"
        script.write_text(
            "#!/bin/sh\n"
            f"while kill -0 {pid} 2>/dev/null; do sleep 0.5; done\n"
            f'rm -rf "{old}"\n'
            f'if mv "{target}" "{old}" && ditto "{new_app}" "{target}"; then\n'
            f'  rm -rf "{old}"\n'
            "else\n"
            f'  rm -rf "{target}"; mv "{old}" "{target}"\n'
            "fi\n"
            f'xattr -dr com.apple.quarantine "{target}" 2>/dev/null\n'
            f'open "{target}"\n'
            f'rm -rf "{stage}"\n', encoding="utf-8")

        def start():
            subprocess.Popen(["/bin/sh", str(script)], start_new_session=True, close_fds=True)
        return start
    raise UpdateError("this installation cannot be updated automatically")


def open_externally(target: str | Path) -> None:
    """Opens a file or URL with the system (Finder, Explorer or the browser)."""
    t = str(target)
    if t.startswith("http"):
        import webbrowser
        webbrowser.open(t)
    elif sys.platform == "darwin":
        subprocess.Popen(["open", t])
    elif sys.platform == "win32":
        os.startfile(t)     # noqa: S606
    else:
        subprocess.Popen(["xdg-open", t])

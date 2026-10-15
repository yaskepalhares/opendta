"""Update check: versions, release parsing, asset choice and the verified download."""

import hashlib
import io
import json
from pathlib import Path

import pytest

from opendta import updater


def test_versions():
    assert updater.is_newer("v0.3.0", "0.2.0")
    assert updater.is_newer("0.2.10", "0.2.9")
    assert not updater.is_newer("v0.2.0", "0.2.0")
    assert not updater.is_newer("0.3.0-beta1", "0.3.0")
    assert updater.is_newer("0.3.0", "0.3.0rc2")
    assert updater.parse_version("1.2") == (1, 2, 0, 1)


RELEASE = {
    "tag_name": "v0.9.0", "body": "## Notes", "html_url": "https://github.com/yaskepalhares/opendta/releases/tag/v0.9.0",
    "assets": [{"name": n, "browser_download_url": "https://example.invalid/" + n, "size": 10,
                "digest": "sha256:" + "0" * 64} for n in updater.ASSETS.values()],
}


class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_fetch_and_assets():
    rel = updater.fetch_latest(opener=lambda req, timeout: _Resp(json.dumps(RELEASE).encode()))
    assert rel.version == "0.9.0" and rel.notes == "## Notes"
    assert updater.asset_for(rel, "windows-installer").name == "OpenDTA-windows-x64-setup.exe"
    assert updater.asset_for(rel, "macos-arm64").name == "OpenDTA-macos-apple-silicon.dmg"
    assert updater.asset_for(rel, "source") is None
    with pytest.raises(updater.UpdateError):
        updater.fetch_latest(opener=lambda req, timeout: (_ for _ in ()).throw(OSError("offline")))


def test_download_checks_size_and_digest(tmp_path: Path):
    payload = b"OpenDTA" * 1000
    src = tmp_path / "src.bin"
    src.write_bytes(payload)
    good = updater.Asset("OpenDTA-test.bin", src.as_uri(), len(payload),
                         "sha256:" + hashlib.sha256(payload).hexdigest())
    seen = []
    out = updater.download(good, tmp_path / "dl", progress=lambda a, b: seen.append(a))
    assert out.read_bytes() == payload and seen[-1] == len(payload)
    bad = updater.Asset("OpenDTA-bad.bin", src.as_uri(), len(payload), "sha256:" + "0" * 64)
    with pytest.raises(updater.UpdateError):
        updater.download(bad, tmp_path / "dl")
    assert not (tmp_path / "dl" / "OpenDTA-bad.bin").exists()
    with pytest.raises(updater.UpdateError):
        updater.download(updater.Asset("x.bin", src.as_uri(), 5), tmp_path / "dl")


def test_source_install_cannot_update_itself(tmp_path):
    assert updater.install_kind() == "source"
    with pytest.raises(updater.UpdateError):
        updater.plan_install("source", tmp_path / "x")


def test_update_dialog_builds(monkeypatch, tmp_path):
    monkeypatch.setenv("OPENDTA_SETTINGS", str(tmp_path / "prefs.ini"))
    QtWidgets = pytest.importorskip("PySide6.QtWidgets")
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])  # noqa: F841
    from opendta.gui.main_window import MainWindow
    from opendta.gui.update_dialog import UpdateDialog
    w = MainWindow()
    d = UpdateDialog(w.updates, updater.parse_release(RELEASE), w)
    assert "0.9.0" in d.findChildren(QtWidgets.QLabel)[0].text()
    w.close()


@pytest.mark.parametrize("fails", [True, False])
def test_check_result_handled_in_gui_thread(monkeypatch, tmp_path, fails):
    """The check runs in a QThread; its result (error message or update dialog) must be shown
    from the GUI thread: a dialog opened in the worker thread aborts OpenDTA on macOS."""
    import threading
    import time

    monkeypatch.setenv("OPENDTA_SETTINGS", str(tmp_path / "prefs.ini"))
    QtWidgets = pytest.importorskip("PySide6.QtWidgets")
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from opendta.gui import update_dialog as UD
    from opendta.gui.main_window import MainWindow

    def fake_fetch():
        if fails:
            raise updater.UpdateError("certificate verify failed")
        return updater.parse_release(RELEASE)

    monkeypatch.setattr(updater, "fetch_latest", fake_fetch)
    seen = []
    monkeypatch.setattr(UD.QMessageBox, "warning", lambda *a, **k: seen.append(threading.current_thread()))
    monkeypatch.setattr(UD.QMessageBox, "information",
                        lambda *a, **k: seen.append(threading.current_thread()))
    monkeypatch.setattr(UD.UpdateDialog, "exec", lambda self: seen.append(threading.current_thread()))
    w = MainWindow()
    w.updates.check_now()
    deadline = time.time() + 10
    while not seen and time.time() < deadline:
        app.processEvents()
        time.sleep(0.01)
    assert seen and seen[0] is threading.main_thread()
    w.updates._thread.wait(5000)
    w.close()


def test_ssl_context_uses_certifi():
    certifi = pytest.importorskip("certifi")
    ctx = updater.ssl_context()
    assert ctx.verify_mode.name == "CERT_REQUIRED"
    assert ctx.cert_store_stats()["x509_ca"] > 50
    assert certifi.where()

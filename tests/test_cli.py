import os
from pathlib import Path

from opendta.cli import run_batch


def test_batch_writes_log(tmp_path: Path):
    do = tmp_path / "prog.do"
    do.write_text('display "ok"\nforvalues i = 1/2 {\n  display `i\'\n}\n', encoding="utf-8")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        rc = run_batch(str(do), [])
    finally:
        os.chdir(cwd)
    assert rc == 0
    log = (tmp_path / "prog.log").read_text(encoding="utf-8")
    assert "log type:  text" in log
    assert '. display "ok"\nok\n' in log
    assert "  2.   display `i'" in log or "  2. display `i'" in log
    assert log.rstrip().endswith("end of do-file")


def test_batch_returns_error_code(tmp_path: Path):
    do = tmp_path / "erro.do"
    do.write_text("foo\n", encoding="utf-8")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        rc = run_batch(str(do), [])
    finally:
        os.chdir(cwd)
    assert rc == 199

import pytest

from opendta.core.output import Capture
from opendta.session import Session


class Runner:
    def __init__(self):
        self.session = Session()
        # the tests check output in Stata's format; OpenDTA's error explanations have their own
        # tests (test_hints.py)
        self.session.settings["hints"] = "off"
        # likewise for the superscript exponent (test_superscript.py)
        self.session.run("quietly set superscript off")
        self.session.run("quietly set numerics stata")
        self.cap = Capture()
        self.session.output.add_listener(self.cap)

    def __call__(self, source: str) -> str:
        """Runs without echo and returns only the output produced."""
        self.cap.clear()
        self.rc = self.session.run(source)
        return self.cap.text

    def eval(self, expr: str):
        return self.session.eval(expr)


@pytest.fixture
def run():
    return Runner()


@pytest.fixture(autouse=True)
def stata_exponents():
    """The tests check Stata's form (e+10); Session() turns the superscript
    exponent on and the state is global, so each test starts with it off."""
    from opendta.core.formats import set_superscript
    from opendta.estimation.numerics import set_precise
    set_superscript(False)
    set_precise(False)
    yield
    set_superscript(False)
    set_precise(False)


@pytest.fixture(autouse=True)
def english_interface():
    """set language is global as well: each test starts and ends in English."""
    from opendta import i18n
    i18n.set_language("en")
    i18n.set_messages("stata")
    yield
    i18n.set_language("en")
    i18n.set_messages("stata")


def pytest_runtest_logreport(report):
    """In CI, each failing test becomes an annotation (the job logs cannot always be read)."""
    import os
    if report.failed and os.environ.get("GITHUB_ACTIONS"):
        msg = str(report.longrepr).splitlines()[-12:]
        print(f"\n::error title=pytest {report.nodeid}::" + "%0A".join(m.replace("%", "%25") for m in msg))

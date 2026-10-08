import pytest

from opendta.core.output import Capture
from opendta.session import Session


class Runner:
    def __init__(self):
        self.session = Session()
        self.cap = Capture()
        self.session.output.add_listener(self.cap)

    def __call__(self, source: str) -> str:
        """Executa sem eco e devolve só a saída produzida."""
        self.cap.clear()
        self.rc = self.session.run(source)
        return self.cap.text

    def eval(self, expr: str):
        return self.session.eval(expr)


@pytest.fixture
def run():
    return Runner()

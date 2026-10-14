"""Ponto de entrada do aplicativo empacotado (PyInstaller).

Sem argumentos abre a interface; com um arquivo (.do ou .dta), abre a
interface e executa ou carrega o arquivo. As opções da linha de comando
(--console, -b) continuam valendo para quem chama o executável de dentro
do pacote pelo Terminal.
"""

import sys

from opendta.cli import main

if __name__ == "__main__":
    # o Finder antigo passa -psn_0_12345 ao abrir o app; não é argumento nosso
    sys.argv = [a for a in sys.argv if not a.startswith("-psn_")]
    sys.exit(main())

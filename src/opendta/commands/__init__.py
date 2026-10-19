"""OpenDTA commands. Importing this package registers all of them."""

from . import adocmd, codebook, putexcel, collapse, combine, data, datatools, delimited, display, egen, estcmds, filecmd, files, graphcmd, helpcmd, matacmd, reshape, summarize, tabulate, infile, inspect, logcmd, matrix, preserve, progcmd, label, macro, program, scalar, syntaxcmd, system, uicmd, nptests, immediate, power, epitab, epiodds, stdize, kappa, roc, rocreg, datautil, datacmds, datasets  # noqa: F401
from .registry import REGISTRY, lookup  # noqa: F401

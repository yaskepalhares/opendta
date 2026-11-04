"""Registration of the model-fitting commands of the Models block (scobit, hetprobit, fracreg,
zip, zinb, tpoisson, tnbreg, gnbreg, binreg...). The logic lives in opendta/estimation; here only
the names and abbreviations are registered."""

from __future__ import annotations

from .registry import lazy_commands

lazy_commands("binmodels", "scobit", "hetprobit", "fracreg", "zip", "zinb", "tpoisson", "tnbreg",
              "gnbreg", "binreg", "glogit")

lazy_commands("censmodels", "intreg", "truncreg", "heckman", "heckprobit", "heckoprobit", "biprobit")

lazy_commands("ivmodels", "ivprobit", "ivtobit", "ivpoisson")

lazy_commands("sysmodels", "constraint", "cnsreg", "eivreg", "vwls", "sureg", "reg3", "rreg")

lazy_commands("choice", "clogit", "rologit", "asclogit", "slogit")

lazy_commands("choiceprobit", "mprobit", "asmprobit", "asroprobit")

lazy_commands("exact", "exlogistic", "expoisson")

lazy_commands("nonlin", "nl", "nlsur")
lazy_commands("gmm", "gmm")

lazy_commands("frontier", "frontier")
lazy_commands("treatment", "etregress", "etpoisson")
lazy_commands("teffects", "teffects", "eteffects", "stteffects")

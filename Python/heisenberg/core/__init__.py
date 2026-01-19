"""
Core module for Heisenberg analysis.

Contains physical constants, the main analysis pipeline, fitting algorithms,
and derived physics calculations.
"""

from heisenberg.core.constants import ASTRO, PHYS, NUMBERS
from heisenberg.core.model import (
    f_zeta,
    f_rpeak,
    f_fluxratiostar,
    f_fluxratiogas,
    f_intfrac,
    gaussian,
)
from heisenberg.core.derived import (
    f_esf,
    f_mdotsf,
    f_mdotfb,
    f_vfb,
    f_vfbr,
    f_etainst,
    f_etaavg,
    f_pzero,
    f_chie,
    f_chier,
    f_chip,
    f_chipr,
)
from heisenberg.core.fitting import fit_kl14, FitResult, pdf_to_values

__all__ = [
    # Constants
    "ASTRO",
    "PHYS",
    "NUMBERS",
    # Model functions
    "f_zeta",
    "f_rpeak",
    "f_fluxratiostar",
    "f_fluxratiogas",
    "f_intfrac",
    "gaussian",
    # Derived quantities
    "f_esf",
    "f_mdotsf",
    "f_mdotfb",
    "f_vfb",
    "f_vfbr",
    "f_etainst",
    "f_etaavg",
    "f_pzero",
    "f_chie",
    "f_chier",
    "f_chip",
    "f_chipr",
    # Fitting
    "fit_kl14",
    "FitResult",
    "pdf_to_values",
]

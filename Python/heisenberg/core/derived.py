"""
Derived Physical Quantities.

This module calculates derived physical quantities from the fitted KL14
model parameters, including star formation efficiencies, feedback velocities,
mass loading factors, and feedback efficiencies.

References:
    Kruijssen & Longmore (2014), MNRAS 439, 3239
    Kruijssen et al. (2018), MNRAS 479, 1866
"""

import numpy as np
from typing import Optional

from heisenberg.core.constants import NUMBERS, ASTRO, PHYS
from heisenberg.core.model import f_zeta, f_rpeak


def f_esf(tgas: float, tdepl: float, fcl: float, fgmc: float) -> float:
    """
    Calculate the star formation efficiency per star formation event.

    Args:
        tgas: Gas phase duration (Myr)
        tdepl: Depletion time (Gyr)
        fcl: Fraction of stellar emission in compact regions
        fgmc: Fraction of gas in giant molecular clouds

    Returns:
        Star formation efficiency (dimensionless, 0-1)
    """
    # Cannot exceed unity
    return min(tgas / (tdepl * 1.0e3), 1.0)


def f_mdotsf(
    tgas: float,
    lambda_: float,
    surfgas: float,
    fgmc: float,
    esf: float
) -> float:
    """
    Calculate the star formation rate per star formation event during gas phase.

    Args:
        tgas: Gas phase duration (Myr)
        lambda_: Mean separation between regions (pc)
        surfgas: Gas surface density (Msun/pc^2)
        fgmc: Fraction of gas in GMCs
        esf: Star formation efficiency

    Returns:
        Star formation rate (Msun/yr)
    """
    # Area of region * gas fraction * efficiency * surface density / time
    mdotsf = np.pi * lambda_**2 * fgmc * esf * surfgas / (4.0 * tgas * 1.0e6)
    return mdotsf


def f_mdotfb(
    tover: float,
    lambda_: float,
    surfgas: float,
    fgmc: float,
    esf: float
) -> float:
    """
    Calculate the mass outflow rate per event during the overlap phase.

    Args:
        tover: Overlap phase duration (Myr)
        lambda_: Mean separation between regions (pc)
        surfgas: Gas surface density (Msun/pc^2)
        fgmc: Fraction of gas in GMCs
        esf: Star formation efficiency

    Returns:
        Mass outflow rate (Msun/yr)
    """
    mdotfb = np.pi * lambda_**2 * fgmc * (1.0 - esf) * surfgas / (4.0 * tover * 1.0e6)
    return max(mdotfb, NUMBERS.tiny)


def f_vfb(tover: float, lambda_: float) -> float:
    """
    Calculate the feedback-driven expansion velocity of ejecta.

    Args:
        tover: Overlap phase duration (Myr)
        lambda_: Mean separation between regions (pc)

    Returns:
        Expansion velocity (km/s)
    """
    # v = 0.5 * lambda / tover, converted to km/s
    vfb = 0.5 * (lambda_ * ASTRO.pc) / (tover * ASTRO.myr) / ASTRO.kms
    return vfb


def f_vfbr(tover: float, rpeakgas: float) -> float:
    """
    Calculate the feedback-driven evaporation velocity of residual gas.

    Uses the peak radius instead of lambda for the characteristic length scale.

    Args:
        tover: Overlap phase duration (Myr)
        rpeakgas: Gas peak radius (pc)

    Returns:
        Evaporation velocity (km/s)
    """
    vfbr = (rpeakgas * ASTRO.pc) / (tover * ASTRO.myr) / ASTRO.kms
    return vfbr


def f_etainst(tgas: float, tover: float, esf: float) -> float:
    """
    Calculate the instantaneous mass loading factor.

    The mass loading factor is the ratio of mass outflow rate to star
    formation rate during the overlap phase.

    Args:
        tgas: Gas phase duration (Myr)
        tover: Overlap phase duration (Myr)
        esf: Star formation efficiency

    Returns:
        Instantaneous mass loading factor (dimensionless)
    """
    etainst = (1.0 - esf) * tgas / (esf * tover)
    return max(etainst, NUMBERS.tiny)


def f_etaavg(esf: float) -> float:
    """
    Calculate the time-integrated mass loading factor.

    Args:
        esf: Star formation efficiency

    Returns:
        Time-integrated mass loading factor (dimensionless)
    """
    etaavg = (1.0 - esf) / esf
    return max(etaavg, NUMBERS.tiny)


def f_pzero(vfbr: float, etaavg: float) -> float:
    """
    Calculate the specific terminal momentum.

    Args:
        vfbr: Feedback evaporation velocity (km/s)
        etaavg: Time-integrated mass loading factor

    Returns:
        Specific terminal momentum (km/s)
    """
    return vfbr * etaavg


def f_chie(tover: float, esf: float, vfb: float, psie: float) -> float:
    """
    Calculate the feedback energy efficiency.

    Args:
        tover: Overlap phase duration (Myr)
        esf: Star formation efficiency
        vfb: Feedback expansion velocity (km/s)
        psie: Energy output rate per unit mass (m^2/s^3)

    Returns:
        Feedback energy efficiency (dimensionless)
    """
    chie = (1.0 - esf) * (vfb * ASTRO.kms)**2 / (2.0 * esf * (tover * ASTRO.myr) * psie)
    return max(chie, NUMBERS.tiny)


def f_chier(tover: float, esf: float, vfbr: float, psie: float) -> float:
    """
    Calculate the feedback energy efficiency using region radius.

    Args:
        tover: Overlap phase duration (Myr)
        esf: Star formation efficiency
        vfbr: Feedback evaporation velocity (km/s)
        psie: Energy output rate per unit mass (m^2/s^3)

    Returns:
        Feedback energy efficiency (dimensionless)
    """
    chier = (1.0 - esf) * (vfbr * ASTRO.kms)**2 / (2.0 * esf * (tover * ASTRO.myr) * psie)
    return max(chier, NUMBERS.tiny)


def f_chip(tover: float, esf: float, vfb: float, psip: float) -> float:
    """
    Calculate the feedback momentum efficiency.

    Args:
        tover: Overlap phase duration (Myr)
        esf: Star formation efficiency
        vfb: Feedback expansion velocity (km/s)
        psip: Momentum output rate per unit mass (m/s^2)

    Returns:
        Feedback momentum efficiency (dimensionless)
    """
    chip = (1.0 - esf) * (vfb * ASTRO.kms) / (esf * (tover * ASTRO.myr) * psip)
    return max(chip, NUMBERS.tiny)


def f_chipr(tover: float, esf: float, vfbr: float, psip: float) -> float:
    """
    Calculate the feedback momentum efficiency using region radius.

    Args:
        tover: Overlap phase duration (Myr)
        esf: Star formation efficiency
        vfbr: Feedback evaporation velocity (km/s)
        psip: Momentum output rate per unit mass (m/s^2)

    Returns:
        Feedback momentum efficiency (dimensionless)
    """
    chipr = (1.0 - esf) * (vfbr * ASTRO.kms) / (esf * (tover * ASTRO.myr) * psip)
    return max(chipr, NUMBERS.tiny)

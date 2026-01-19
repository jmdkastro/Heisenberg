"""
KL14 Uncertainty Principle Model Functions.

This module implements the core mathematical model from Kruijssen & Longmore (2014)
for constraining the timescales of star formation and feedback in galaxies.

The key parameters fitted are:
- tgas: Duration of the gas phase (Myr)
- tover: Duration of the overlap phase (Myr)
- lambda: Mean separation between independent star-forming regions (pc)

References:
    Kruijssen & Longmore (2014), MNRAS 439, 3239
    Kruijssen et al. (2018), MNRAS 479, 1866
"""

import numpy as np
from typing import Tuple, Optional
from scipy.special import erf

from heisenberg.core.constants import NUMBERS


def f_zeta(surfcontrast: float, ttot: float, tref: float, peak_prof: int) -> float:
    """
    Calculate the tracer peak concentration parameter zeta.

    Args:
        surfcontrast: Surface brightness contrast of peaks
        ttot: Total timeline duration (tgas + tstar - tover) in Myr
        tref: Reference timescale (tgas or tstar) in Myr
        peak_prof: Peak profile model (0=point, 1=disc, 2=Gaussian)

    Returns:
        Concentration parameter zeta (dimensionless)
    """
    surfcontrast = max(surfcontrast, 1.0)

    if peak_prof == 0:
        return NUMBERS.tiny
    elif peak_prof == 1:
        return np.sqrt((ttot / tref) / (surfcontrast * (1.0 + ttot / tref) - 1.0))
    elif peak_prof == 2:
        return np.sqrt((ttot / tref) / (2.0 * surfcontrast * (1.0 + ttot / tref) - 2.0))
    else:
        raise ValueError(f"Invalid peak_prof: {peak_prof}. Must be 0, 1, or 2.")


def f_rpeak(zeta: float, lambda_: float) -> float:
    """
    Calculate the tracer peak radius.

    Args:
        zeta: Concentration parameter from f_zeta
        lambda_: Mean separation between regions (pc)

    Returns:
        Peak radius in pc
    """
    return 0.5 * zeta * lambda_


def f_fluxratiostar(
    tgas: float,
    tstar: float,
    tover: float,
    laps: np.ndarray,
    lambda_: float,
    beta_gas: float,
    surfcontrasts: float,
    surfcontrastg: float,
    peak_prof: int
) -> np.ndarray:
    """
    Calculate the stellar flux ratio (depletion time ratio) at each aperture.

    This is the ratio of depletion time measured around stellar peaks to the
    galactic average depletion time.

    Args:
        tgas: Gas phase duration (Myr)
        tstar: Stellar tracer visibility time (Myr)
        tover: Overlap phase duration (Myr)
        laps: Array of aperture sizes (pc)
        lambda_: Mean separation between regions (pc)
        beta_gas: Gas flux ratio during overlap
        surfcontrasts: Stellar surface brightness contrast
        surfcontrastg: Gas surface brightness contrast
        peak_prof: Peak profile model (0=point, 1=disc, 2=Gaussian)

    Returns:
        Array of flux ratios at each aperture size
    """
    ttotal = tgas + tstar - tover
    zetas = f_zeta(surfcontrasts, ttotal, tstar, peak_prof)
    zetag = f_zeta(surfcontrastg, ttotal, tgas, peak_prof)
    rpeaks = f_rpeak(zetas, lambda_)
    rpeakg = f_rpeak(zetag, lambda_)

    laps = np.atleast_1d(laps)
    nl = len(laps)
    surfratios = np.zeros(nl)
    surfratiog = np.zeros(nl)

    for i in range(nl):
        if peak_prof <= 1:
            # Point or disc model
            if rpeaks < 0.5 * laps[i]:
                surfratios[i] = ttotal / tstar * (lambda_ / laps[i])**2
            else:
                surfratios[i] = surfcontrasts
            if rpeakg < 0.5 * laps[i]:
                surfratiog[i] = ttotal / tgas * (lambda_ / laps[i])**2
            else:
                surfratiog[i] = surfcontrastg
        elif peak_prof == 2:
            # Gaussian model
            surfratios[i] = ttotal / tstar * (lambda_ / laps[i])**2 * (
                1.0 - np.exp(-laps[i]**2 / 8.0 / rpeaks**2)
            )
            surfratiog[i] = ttotal / tgas * (lambda_ / laps[i])**2 * (
                1.0 - np.exp(-laps[i]**2 / 8.0 / rpeakg**2)
            )

    others = 1.0
    gas = (beta_gas * tover / tstar / (1.0 + (beta_gas - 1.0) * tover / tgas)) * surfratiog
    star = surfratios
    tdeplstar = (gas + others) / (star + others)

    return tdeplstar


def f_fluxratiogas(
    tgas: float,
    tstar: float,
    tover: float,
    lapg: np.ndarray,
    lambda_: float,
    beta_star: float,
    surfcontrasts: float,
    surfcontrastg: float,
    peak_prof: int
) -> np.ndarray:
    """
    Calculate the gas flux ratio (depletion time ratio) at each aperture.

    This is the ratio of depletion time measured around gas peaks to the
    galactic average depletion time.

    Args:
        tgas: Gas phase duration (Myr)
        tstar: Stellar tracer visibility time (Myr)
        tover: Overlap phase duration (Myr)
        lapg: Array of aperture sizes (pc)
        lambda_: Mean separation between regions (pc)
        beta_star: Stellar flux ratio during overlap
        surfcontrasts: Stellar surface brightness contrast
        surfcontrastg: Gas surface brightness contrast
        peak_prof: Peak profile model (0=point, 1=disc, 2=Gaussian)

    Returns:
        Array of flux ratios at each aperture size
    """
    ttotal = tgas + tstar - tover
    zetas = f_zeta(surfcontrasts, ttotal, tstar, peak_prof)
    zetag = f_zeta(surfcontrastg, ttotal, tgas, peak_prof)
    rpeaks = f_rpeak(zetas, lambda_)
    rpeakg = f_rpeak(zetag, lambda_)

    lapg = np.atleast_1d(lapg)
    nl = len(lapg)
    surfratios = np.zeros(nl)
    surfratiog = np.zeros(nl)

    for i in range(nl):
        if peak_prof <= 1:
            # Point or disc model
            if rpeaks < 0.5 * lapg[i]:
                surfratios[i] = ttotal / tstar * (lambda_ / lapg[i])**2
            else:
                surfratios[i] = surfcontrasts
            if rpeakg < 0.5 * lapg[i]:
                surfratiog[i] = ttotal / tgas * (lambda_ / lapg[i])**2
            else:
                surfratiog[i] = surfcontrastg
        elif peak_prof == 2:
            # Gaussian model
            surfratios[i] = ttotal / tstar * (lambda_ / lapg[i])**2 * (
                1.0 - np.exp(-lapg[i]**2 / 8.0 / rpeaks**2)
            )
            surfratiog[i] = ttotal / tgas * (lambda_ / lapg[i])**2 * (
                1.0 - np.exp(-lapg[i]**2 / 8.0 / rpeakg**2)
            )

    others = 1.0
    gas = surfratiog
    star = (beta_star * tover / tgas / (1.0 + (beta_star - 1.0) * tover / tstar)) * surfratios
    tdeplgas = (gas + others) / (star + others)

    return tdeplgas


def f_intfrac(var: np.ndarray, r1: float, r2: float) -> float:
    """
    Calculate the integrated fraction of a Gaussian profile between two radii.

    Used for calculating the flux enclosed within apertures for Gaussian
    peak profiles.

    Args:
        var: Array of Gaussian parameters [amplitude, center, sigma]
        r1: Inner radius
        r2: Outer radius

    Returns:
        Fraction of total flux between r1 and r2
    """
    a0 = var[0]  # Gaussian normalisation
    a1 = var[1]  # Gaussian x offset
    a2 = var[2]  # Gaussian dispersion

    # Surface density at inner and outer edges
    sigma1 = a0 * np.exp(-0.5 * ((r1 - a1) / a2)**2)
    sigma2 = a0 * np.exp(-0.5 * ((r2 - a1) / a2)**2)

    # Error functions at inner and outer edges
    erf1 = erf((r1 - a1) / (np.sqrt(2.0) * a2))
    erf2 = erf((r2 - a1) / (np.sqrt(2.0) * a2))

    # Integral of surface density between inner and outer edge
    int1 = a2**2 * (sigma1 - sigma2) + np.sqrt(np.pi / 2.0) * a0 * a1 * a2 * (erf2 - erf1)
    # Integral of surface density between inner edge and infinity
    int2 = a2**2 * sigma1 + np.sqrt(np.pi / 2.0) * a0 * a1 * a2 * (1.0 - erf1)

    # Fraction covered within outer edge
    return int1 / int2


def gaussian(x: np.ndarray, params: np.ndarray) -> np.ndarray:
    """
    Evaluate a Gaussian function.

    Args:
        x: Input array of x values
        params: Array of [amplitude, center, sigma]

    Returns:
        Gaussian values at each x
    """
    a0, a1, a2 = params
    return a0 * np.exp(-0.5 * ((x - a1) / a2)**2)

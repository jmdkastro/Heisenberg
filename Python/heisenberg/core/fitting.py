"""
KL14 Model Fitting.

This module implements the chi-squared fitting procedure for the KL14
uncertainty principle model. It performs a 3D grid search over:
- tgas: Gas phase duration
- tover: Overlap phase duration
- lambda: Mean separation between regions

The fitting is iteratively refined by zooming in on the best-fit region.

References:
    Kruijssen & Longmore (2014), MNRAS 439, 3239
    Kruijssen et al. (2018), MNRAS 479, 1866
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional, List
from scipy.special import erf

from heisenberg.core.constants import NUMBERS
from heisenberg.core.model import f_fluxratiostar, f_fluxratiogas


@dataclass
class FitResult:
    """Results from the KL14 model fitting procedure."""

    # Best-fit values
    tgas: float  # Gas phase duration (Myr)
    tover: float  # Overlap phase duration (Myr)
    lambda_: float  # Mean separation (pc)
    tstar: float  # Stellar tracer duration (Myr)

    # Uncertainties (1-sigma)
    tgas_errmin: float
    tgas_errmax: float
    tover_errmin: float
    tover_errmax: float
    lambda_errmin: float
    lambda_errmax: float

    # Fit quality
    chi2_min: float  # Minimum reduced chi-squared
    ndof: int  # Degrees of freedom

    # Best-fit auxiliary parameters
    beta_star: float
    beta_gas: float
    surfcon_star: float
    surfcon_gas: float

    # Probability arrays (for plotting)
    tgas_arr: np.ndarray
    tover_arr: np.ndarray
    lambda_arr: np.ndarray
    prob_tgas: np.ndarray
    prob_tover: np.ndarray
    prob_lambda: np.ndarray

    # Correlations
    corr_tgas_tover: float
    corr_tgas_lambda: float
    corr_tover_lambda: float


def pdf_to_values(
    variable: np.ndarray,
    dvariable: np.ndarray,
    cumpdf: np.ndarray,
    refvalue: float
) -> Tuple[float, float, float]:
    """
    Convert a PDF to best-fit value with asymmetric error bars.

    The errors are set by the 1-sigma percentiles relative to the
    reference (best-fit) value.

    Args:
        variable: Array of parameter values
        dvariable: Array of bin widths
        cumpdf: Cumulative probability distribution
        refvalue: Reference (best-fit) value

    Returns:
        Tuple of (expectation value, errmin, errmax)
    """
    # Gaussian 1-sigma percentiles
    errminfrac_gaussian = 1.0 - 0.5 * (1.0 + erf(1.0 / np.sqrt(2.0)))
    errmaxfrac_gaussian = 1.0 - errminfrac_gaussian

    # Find probability at reference value
    refprob = np.interp(refvalue, variable, cumpdf)
    refprob = np.clip(refprob, NUMBERS.tiny, 1.0 - NUMBERS.tiny)
    if not np.isfinite(refprob):
        refprob = 0.5

    # Calculate asymmetric error fractions
    errminfrac = 2.0 * errminfrac_gaussian * refprob
    errmaxfrac = 2.0 * (1.0 - errmaxfrac_gaussian) * refprob + 2.0 * (errmaxfrac_gaussian - 0.5)

    # Calculate PDF from cumulative
    nvar = len(variable)
    pdf = np.zeros(nvar)
    pdf[0] = cumpdf[0] / dvariable[0]
    for i in range(1, nvar):
        pdf[i] = (cumpdf[i] - cumpdf[i - 1]) / dvariable[i]

    # Expectation value
    evalue = np.sum(variable * pdf * dvariable)

    # Percentiles
    minperc = np.interp(errminfrac, cumpdf, variable)
    maxperc = np.interp(errmaxfrac, cumpdf, variable)

    if minperc != maxperc:
        errmin = max(0.0, refvalue - minperc)
        errmax = max(0.0, maxperc - refvalue)
    else:
        errmin = NUMBERS.tiny
        errmax = NUMBERS.tiny

    return evalue, errmin, errmax


def fit_kl14(
    fluxratio_star: np.ndarray,
    fluxratio_gas: np.ndarray,
    err_star_log: np.ndarray,
    err_gas_log: np.ndarray,
    tstariso: float,
    beta_star: np.ndarray,
    beta_gas: np.ndarray,
    fstarover: np.ndarray,
    fgasover: np.ndarray,
    apertures_star: np.ndarray,
    apertures_gas: np.ndarray,
    surfcontrasts: np.ndarray,
    surfcontrastg: np.ndarray,
    peak_prof: int,
    tstar_incl: bool,
    tgasmini: float,
    tgasmaxi: float,
    tovermini: float,
    ndepth: int = 4,
    ntry: int = 101,
    nfitstar: Optional[np.ndarray] = None,
    nfitgas: Optional[np.ndarray] = None,
    verbose: bool = True
) -> FitResult:
    """
    Fit the KL14 uncertainty principle model to observed flux ratios.

    This performs a 3D chi-squared minimization over tgas, tover, and lambda
    using an iterative grid refinement approach.

    Args:
        fluxratio_star: Observed stellar flux ratios at each aperture
        fluxratio_gas: Observed gas flux ratios at each aperture
        err_star_log: Logarithmic errors on stellar flux ratios
        err_gas_log: Logarithmic errors on gas flux ratios
        tstariso: Reference stellar tracer timescale (Myr)
        beta_star: Stellar flux ratio during overlap vs overlap fraction
        beta_gas: Gas flux ratio during overlap vs overlap fraction
        fstarover: Overlap fraction array for beta_star interpolation
        fgasover: Overlap fraction array for beta_gas interpolation
        apertures_star: Aperture sizes for stellar measurements (pc)
        apertures_gas: Aperture sizes for gas measurements (pc)
        surfcontrasts: Stellar surface brightness contrast at each aperture
        surfcontrastg: Gas surface brightness contrast at each aperture
        peak_prof: Peak profile model (0=point, 1=disc, 2=Gaussian)
        tstar_incl: If True, tstariso includes overlap phase
        tgasmini: Minimum tgas for fitting (Myr)
        tgasmaxi: Maximum tgas for fitting (Myr)
        tovermini: Minimum tover for fitting (Myr)
        ndepth: Maximum number of refinement iterations
        ntry: Number of grid points per dimension
        nfitstar: Number of peaks at each aperture (stellar), for weighting
        nfitgas: Number of peaks at each aperture (gas), for weighting
        verbose: Print progress messages

    Returns:
        FitResult containing best-fit parameters and uncertainties
    """
    napertures = len(apertures_star)
    if len(apertures_gas) != napertures:
        raise ValueError("apertures_star and apertures_gas must have same length")

    # Find valid apertures (finite flux ratios)
    use = np.where(
        np.isfinite(np.log10(fluxratio_star[:napertures])) &
        np.isfinite(np.log10(fluxratio_gas[:napertures]))
    )[0]

    nvar = 3  # tgas, tover, lambda

    # Set up weights based on number of peaks
    if nfitstar is None:
        nfitstar = np.ones(napertures)
    if nfitgas is None:
        nfitgas = np.ones(napertures)

    nfit = np.sum(nfitstar[use]) + np.sum(nfitgas[use])
    ndeg = nfit - nvar
    weights = np.concatenate([nfitstar[use], nfitgas[use]])
    weights = weights / np.mean(weights)

    # Set initial fitting ranges
    lambdamini = 0.3 * min(np.min(apertures_star[use]), np.min(apertures_gas[use]))
    lambdamaxi = 3.0 * max(np.max(apertures_star[use]), np.max(apertures_gas[use]))

    if tstar_incl:
        tovermaxi = min(tstariso, tgasmaxi)
    else:
        tovermaxi = tgasmaxi

    tgasmin, tgasmax = tgasmini, tgasmaxi
    tovermin, tovermax = tovermini, tovermaxi
    lambdamin, lambdamax = lambdamini, lambdamaxi

    # Iterative grid refinement
    for n in range(ndepth):
        if verbose:
            print(f"     ==> fitting iteration {n + 1}/{ndepth}")

        # Create logarithmically-spaced parameter arrays
        tgasarr = tgasmin * (tgasmax / tgasmin) ** ((np.arange(ntry) + 0.5) / ntry)
        toverarr = tovermin * (tovermax / tovermin) ** ((np.arange(ntry) + 0.5) / ntry)
        lambdaarr = lambdamin * (lambdamax / lambdamin) ** ((np.arange(ntry) + 0.5) / ntry)

        # Bin widths for probability normalization
        dtgas = tgasarr * ((tgasmax / tgasmin) ** (1.0 / (2.0 * ntry)) -
                           (tgasmax / tgasmin) ** (-1.0 / (2.0 * ntry)))
        dtover = toverarr * ((tovermax / tovermin) ** (1.0 / (2.0 * ntry)) -
                             (tovermax / tovermin) ** (-1.0 / (2.0 * ntry)))
        dlambda = lambdaarr * ((lambdamax / lambdamin) ** (1.0 / (2.0 * ntry)) -
                               (lambdamax / lambdamin) ** (-1.0 / (2.0 * ntry)))

        # Interpolate surface contrasts to lambda values
        surfs_use = np.maximum(
            1.0,
            10.0 ** np.interp(np.log10(lambdaarr), np.log10(apertures_star), np.log10(surfcontrasts))
        )
        surfg_use = np.maximum(
            1.0,
            10.0 ** np.interp(np.log10(lambdaarr), np.log10(apertures_gas), np.log10(surfcontrastg))
        )

        # Initialize chi-squared array
        chi2 = np.full((ntry, ntry, ntry), NUMBERS.huge)

        # Grid search
        for i in range(ntry):
            for j in range(ntry):
                # Get tstar for this (tgas, tover) combination
                if tstar_incl:
                    tstaruse = tstariso
                else:
                    tstaruse = tstariso + toverarr[j]

                # Interpolate beta values
                beta_star_use = np.interp(toverarr[j] / tstaruse, fstarover, beta_star)
                beta_gas_use = np.interp(toverarr[j] / tgasarr[i], fgasover, beta_gas)

                for k in range(ntry):
                    # Check physical constraint: tover <= tgas
                    if toverarr[j] <= tgasarr[i]:
                        # Calculate model flux ratios
                        fr_star_th = f_fluxratiostar(
                            tgasarr[i], tstaruse, toverarr[j],
                            apertures_star, lambdaarr[k], beta_gas_use,
                            surfs_use[k], surfg_use[k], peak_prof
                        )
                        fr_gas_th = f_fluxratiogas(
                            tgasarr[i], tstaruse, toverarr[j],
                            apertures_gas, lambdaarr[k], beta_star_use,
                            surfs_use[k], surfg_use[k], peak_prof
                        )

                        # Calculate chi-squared
                        diff_star = (np.log10(fluxratio_star[use] / fr_star_th[use]))**2 / err_star_log[use]**2
                        diff_gas = (np.log10(fluxratio_gas[use] / fr_gas_th[use]))**2 / err_gas_log[use]**2
                        diff = np.concatenate([diff_star, diff_gas])
                        chi2[i, j, k] = np.sum(diff * weights)

        # Find minimum chi-squared
        redchi2 = chi2 / ndeg
        minchi2 = np.min(redchi2)
        idx = np.unravel_index(np.argmin(redchi2), redchi2.shape)
        ibest, jbest, kbest = idx

        # Calculate probability distribution
        prob = np.exp(-chi2 / 2.0)

        # Create 3D arrays for bin widths
        dtgasarr = np.zeros((ntry, ntry, ntry))
        dtoverarr = np.zeros((ntry, ntry, ntry))
        dlambdaarr = np.zeros((ntry, ntry, ntry))
        for i in range(ntry):
            dtgasarr[i, :, :] = dtgas[i]
            dtoverarr[:, i, :] = dtover[i]
            dlambdaarr[:, :, i] = dlambda[i]

        # Normalize probability
        probcst = 1.0 / np.sum(prob * dtgasarr * dtoverarr * dlambdaarr)
        probnorm = prob * probcst

        # Marginalize to get 1D PDFs
        probtgas = np.sum(probnorm * dtoverarr * dlambdaarr, axis=(1, 2))
        probtover = np.sum(probnorm * dtgasarr * dlambdaarr, axis=(0, 2))
        problambda = np.sum(probnorm * dtgasarr * dtoverarr, axis=(0, 1))

        # Calculate means and correlations
        meantgas = np.sum(tgasarr * probtgas * dtgas)
        meantover = np.sum(toverarr * probtover * dtover)
        meanlambda = np.sum(lambdaarr * problambda * dlambda)

        stdevtgas = np.sqrt(np.sum((tgasarr - meantgas)**2 * probtgas * dtgas))
        stdevtover = np.sqrt(np.sum((toverarr - meantover)**2 * probtover * dtover))
        stdevlambda = np.sqrt(np.sum((lambdaarr - meanlambda)**2 * problambda * dlambda))

        # Create 3D coordinate arrays for correlation calculation
        tgasarr3d = np.zeros((ntry, ntry, ntry))
        toverarr3d = np.zeros((ntry, ntry, ntry))
        lambdaarr3d = np.zeros((ntry, ntry, ntry))
        for i in range(ntry):
            tgasarr3d[i, :, :] = tgasarr[i]
            toverarr3d[:, i, :] = toverarr[i]
            lambdaarr3d[:, :, i] = lambdaarr[i]

        meantgastover = np.sum(tgasarr3d * toverarr3d * probnorm * dtgasarr * dtoverarr * dlambdaarr)
        meantgaslambda = np.sum(tgasarr3d * lambdaarr3d * probnorm * dtgasarr * dtoverarr * dlambdaarr)
        meantoverlambda = np.sum(toverarr3d * lambdaarr3d * probnorm * dtgasarr * dtoverarr * dlambdaarr)

        corrtgastover = (meantgastover - meantgas * meantover) / (stdevtgas * stdevtover) if stdevtgas * stdevtover > 0 else 0
        corrtgaslambda = (meantgaslambda - meantgas * meanlambda) / (stdevtgas * stdevlambda) if stdevtgas * stdevlambda > 0 else 0
        corrtoverlambda = (meantoverlambda - meantover * meanlambda) / (stdevtover * stdevlambda) if stdevtover * stdevlambda > 0 else 0

        # Best-fit values
        tgas = tgasarr[ibest]
        tover = toverarr[jbest]
        lambda_ = lambdaarr[kbest]

        if tstar_incl:
            tstar = tstariso
        else:
            tstar = tstariso + tover

        # Best-fit auxiliary parameters
        beta_star_best = np.interp(tover / tstar, fstarover, beta_star)
        beta_gas_best = np.interp(tover / tgas, fgasover, beta_gas)
        surfs_best = 10.0 ** np.interp(np.log10(lambda_), np.log10(apertures_star), np.log10(surfcontrasts))
        surfg_best = 10.0 ** np.interp(np.log10(lambda_), np.log10(apertures_gas), np.log10(surfcontrastg))

        # Calculate uncertainties from cumulative PDFs
        probtgascum = np.cumsum(probtgas * dtgas)
        probtovercum = np.cumsum(probtover * dtover)
        problambdacum = np.cumsum(problambda * dlambda)

        _, tgas_errmin, tgas_errmax = pdf_to_values(tgasarr, dtgas, probtgascum, tgas)
        _, tover_errmin, tover_errmax = pdf_to_values(toverarr, dtover, probtovercum, tover)
        _, lambda_errmin, lambda_errmax = pdf_to_values(lambdaarr, dlambda, problambdacum, lambda_)

        # Check if we should stop refining
        edge = 4
        tol = 0.01

        # Check if all 1-sigma regions are well within the grid
        if (ibest >= edge and ibest <= ntry - edge - 1 and
            jbest >= edge and jbest <= ntry - edge - 1 and
            kbest >= edge and kbest <= ntry - edge - 1):

            # Check if fitting range has converged
            tgasminold, tgasmaxold = tgasmin, tgasmax
            toverminold, tovermaxold = tovermin, tovermax
            lambdaminold, lambdamaxold = lambdamin, lambdamax

            # Zoom in on the chi-squared minimum region
            deltachi = 14.2 / ndeg  # 3-sigma for 3 parameters

            iredchi = np.min(redchi2, axis=(1, 2))
            jredchi = np.min(redchi2, axis=(0, 2))
            kredchi = np.min(redchi2, axis=(0, 1))

            imin = max(0, np.min(np.where(iredchi <= minchi2 + deltachi)[0]) - 1)
            imax = min(ntry - 1, np.max(np.where(iredchi <= minchi2 + deltachi)[0]) + 1)
            jmin = max(0, np.min(np.where(jredchi <= minchi2 + deltachi)[0]) - 1)
            jmax = min(ntry - 1, np.max(np.where(jredchi <= minchi2 + deltachi)[0]) + 1)
            kmin = max(0, np.min(np.where(kredchi <= minchi2 + deltachi)[0]) - 1)
            kmax = min(ntry - 1, np.max(np.where(kredchi <= minchi2 + deltachi)[0]) + 1)

            tgasmin = tgasarr[imin]
            tgasmax = tgasarr[imax]
            tovermin = min(toverarr[jmin], tgasmin)
            tovermax = toverarr[jmax]
            lambdamin = lambdaarr[kmin]
            lambdamax = lambdaarr[kmax]

            # Check for convergence
            if (abs(np.log10(tgasminold / tgasmin)) <= tol and
                abs(np.log10(tgasmaxold / tgasmax)) <= tol and
                abs(np.log10(toverminold / tovermin)) <= tol and
                abs(np.log10(tovermaxold / tovermax)) <= tol and
                abs(np.log10(lambdaminold / lambdamin)) <= tol and
                abs(np.log10(lambdamaxold / lambdamax)) <= tol):
                if verbose:
                    print("     ==> fitting converged")
                break
        else:
            # Avoid edge traps by extending the range
            if ibest < edge:
                tgasmin = max(tgasmin**2 / tgasmax, tgasmini)
            if ibest > ntry - edge - 1:
                tgasmax = min(tgasmax**2 / tgasmin, tgasmaxi)
            if jbest < edge:
                tovermin = max(2 * tovermin - tovermax, tovermini)
            if jbest > ntry - edge - 1:
                tovermax = min(2 * tovermax - tovermin, tovermaxi)
            if kbest < edge:
                lambdamin = max(lambdamin**2 / lambdamax, lambdamini)
            if kbest > ntry - edge - 1:
                lambdamax = min(lambdamax**2 / lambdamin, lambdamaxi)

    return FitResult(
        tgas=tgas,
        tover=tover,
        lambda_=lambda_,
        tstar=tstar,
        tgas_errmin=tgas_errmin,
        tgas_errmax=tgas_errmax,
        tover_errmin=tover_errmin,
        tover_errmax=tover_errmax,
        lambda_errmin=lambda_errmin,
        lambda_errmax=lambda_errmax,
        chi2_min=minchi2,
        ndof=int(ndeg),
        beta_star=beta_star_best,
        beta_gas=beta_gas_best,
        surfcon_star=surfs_best,
        surfcon_gas=surfg_best,
        tgas_arr=tgasarr,
        tover_arr=toverarr,
        lambda_arr=lambdaarr,
        prob_tgas=probtgas,
        prob_tover=probtover,
        prob_lambda=problambda,
        corr_tgas_tover=corrtgastover,
        corr_tgas_lambda=corrtgaslambda,
        corr_tover_lambda=corrtoverlambda,
    )

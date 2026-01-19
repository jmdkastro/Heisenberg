"""
Fourier Filtering Corrections.

Implements corrections for flux loss due to filtering:
- Flux loss correction (qcon): For isolated Gaussian sources
- Overlap correction (qover): For overlapping peaks

Based on Hygate et al. (2019) Equations 31 and related calibrations.

References:
    Hygate et al. (2019)
    IDL: fourier_flux_loss_correction.pro, fourier_overlap_correction.pro
"""

import numpy as np
from typing import Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum

from .filters import FilterKernel


# ============================================================================
# Sigmoid Functions
# ============================================================================

def symmetric_sigmoidal(
    x: Union[float, np.ndarray],
    a: float,
    b: float,
    c: float,
    d: float = 1.0,
) -> Union[float, np.ndarray]:
    """
    Symmetric sigmoidal function.

    Formula:
        f(x) = d + (a - d) / (1 + (x/c)^b)

    Args:
        x: Input value(s)
        a: Asymptotic minimum
        b: Slope parameter
        c: Inflection point
        d: Asymptotic maximum (default 1.0 for unity sigmoid)

    Returns:
        Sigmoidal output
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        result = d + (a - d) / (1.0 + (x / c)**b)
        # Handle x=0 case
        if np.isscalar(x):
            if x == 0:
                result = a
        else:
            result = np.where(x == 0, a, result)
    return result


def asymmetric_sigmoidal(
    x: Union[float, np.ndarray],
    a: float,
    b: float,
    c: float,
    d: float,
) -> Union[float, np.ndarray]:
    """
    Asymmetric sigmoidal function (same formula as symmetric but with general d).

    Formula:
        f(x) = d + (a - d) / (1 + (x/c)^b)

    Args:
        x: Input value(s)
        a: Lower asymptote
        b: Slope parameter
        c: Inflection point
        d: Upper asymptote

    Returns:
        Sigmoidal output
    """
    return symmetric_sigmoidal(x, a, b, c, d)


# ============================================================================
# Calibrated Parameters
# ============================================================================

# Flux loss correction parameters (Gaussian kernel)
# From Hygate et al. (2019) calibration
FLUX_LOSS_PARAMS_GAUSSIAN = {
    'a': -0.0159072,
    'b': 1.68958,
    'c': 4.85505,
}

# Overlap correction parameters (Gaussian kernel)
OVERLAP_PARAMS_GAUSSIAN = {
    'a': -0.35443467,
    'b': 3.1018254,
    'c': 1.40684605,
    'd': 0.98407271,
}

# Covariance matrix for overlap correction error propagation
# Rows/columns correspond to [a, b, c, d]
OVERLAP_COVARIANCE_GAUSSIAN = np.array([
    [1.80773315e-02, -2.82773696e-02, -4.61227583e-03, -5.61614628e-04],
    [-2.82773696e-02, 1.24531882e-01, 2.35316254e-02, 3.31696538e-04],
    [-4.61227583e-03, 2.35316254e-02, 8.49772673e-03, -7.35095656e-04],
    [-5.61614628e-04, 3.31696538e-04, -7.35095656e-04, 3.01127749e-04],
])


# ============================================================================
# Flux Loss Correction
# ============================================================================

def flux_loss_correction(
    cut_ratio: Union[float, np.ndarray],
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
) -> Union[float, np.ndarray]:
    """
    Calculate flux loss correction factor for isolated sources.

    When filtering removes flux from compact sources, this correction
    factor (qcon) compensates for the loss. Applied as:
        corrected_flux = measured_flux / qcon

    Args:
        cut_ratio: Ratio of filter cut length to source size
            (cut_length / peak_radius)
        kernel: Filter kernel type (only Gaussian implemented)

    Returns:
        Correction factor qcon (typically < 1)

    Notes:
        - Based on Hygate et al. (2019) Equation 31
        - Uses unity symmetric sigmoidal (tends to 1 as cut_ratio increases)
        - For Butterworth and Ideal kernels, returns 1.0 (not yet calibrated)
    """
    if isinstance(kernel, str):
        kernel = FilterKernel(kernel.lower())

    if kernel == FilterKernel.GAUSSIAN:
        params = FLUX_LOSS_PARAMS_GAUSSIAN
        qcon = symmetric_sigmoidal(
            cut_ratio,
            a=params['a'],
            b=params['b'],
            c=params['c'],
            d=1.0,  # Unity sigmoid
        )
    else:
        # Not yet calibrated for other kernels
        qcon = 1.0 if np.isscalar(cut_ratio) else np.ones_like(cut_ratio)

    return qcon


# ============================================================================
# Overlap Correction
# ============================================================================

def overlap_correction(
    dist_stat: Union[float, np.ndarray],
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
) -> Union[float, np.ndarray]:
    """
    Calculate overlap correction factor for crowded peaks.

    When peaks overlap, filtering causes additional flux loss.
    This correction factor (qover) compensates. Applied as:
        corrected_flux = measured_flux / qover

    Args:
        dist_stat: Distance statistic (median nearest-neighbor distance
            normalized by filter cut length)
        kernel: Filter kernel type (only Gaussian implemented)

    Returns:
        Correction factor qover (clamped to max 1.0)

    Notes:
        - Uses asymmetric sigmoidal function
        - Returns 1.0 for Butterworth and Ideal (not yet calibrated)
    """
    if isinstance(kernel, str):
        kernel = FilterKernel(kernel.lower())

    if kernel == FilterKernel.GAUSSIAN:
        params = OVERLAP_PARAMS_GAUSSIAN
        qover = asymmetric_sigmoidal(
            dist_stat,
            a=params['a'],
            b=params['b'],
            c=params['c'],
            d=params['d'],
        )
        # Clamp to maximum of 1.0
        if np.isscalar(qover):
            qover = min(qover, 1.0)
        else:
            qover = np.minimum(qover, 1.0)
    else:
        # Not yet calibrated for other kernels
        qover = 1.0 if np.isscalar(dist_stat) else np.ones_like(dist_stat)

    return qover


def overlap_correction_sigma(
    dist_stat: float,
    dist_stat_sigma: float = 0.0,
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
) -> float:
    """
    Calculate 1-sigma uncertainty on overlap correction factor.

    Uses error propagation through the sigmoidal function with
    the calibrated covariance matrix.

    Args:
        dist_stat: Distance statistic value
        dist_stat_sigma: Uncertainty on distance statistic
        kernel: Filter kernel type

    Returns:
        1-sigma uncertainty on qover

    Notes:
        Based on Hughes & Hase (2010) error propagation method.
    """
    if isinstance(kernel, str):
        kernel = FilterKernel(kernel.lower())

    if kernel != FilterKernel.GAUSSIAN:
        return 0.0

    params = OVERLAP_PARAMS_GAUSSIAN
    a, b, c, d = params['a'], params['b'], params['c'], params['d']
    x = dist_stat
    cov = OVERLAP_COVARIANCE_GAUSSIAN

    # Avoid division by zero
    if x == 0 or c == 0:
        return 0.0

    # Compute ratio and denominator terms
    ratio = x / c
    ratio_b = ratio**b
    denom = 1.0 + ratio_b
    denom_sq = denom**2

    # Partial derivatives
    df_da = 1.0 / denom
    df_db = -((a - d) * ratio_b * np.log(ratio)) / denom_sq if ratio > 0 else 0.0
    df_dc = (b * (a - d) * ratio_b) / (c * denom_sq)
    df_dd = 1.0 - 1.0 / denom
    df_dx = (b * (d - a) * ratio_b) / (x * denom_sq)

    # Build gradient vector (5 elements: a, b, c, d, x)
    partials = np.array([df_da, df_db, df_dc, df_dd, df_dx])

    # Build 5x5 covariance matrix (add x uncertainty on diagonal)
    cov_full = np.zeros((5, 5))
    cov_full[:4, :4] = cov
    cov_full[4, 4] = dist_stat_sigma**2

    # Error propagation: sigma² = sum_ij (∂f/∂pi * ∂f/∂pj * cov_ij)
    variance = 0.0
    for i in range(5):
        for j in range(5):
            variance += partials[i] * partials[j] * cov_full[i, j]

    return np.sqrt(max(variance, 0.0))


# ============================================================================
# Combined Correction
# ============================================================================

@dataclass
class CorrectionFactors:
    """
    Combined correction factors for Fourier filtering.

    Attributes:
        qcon: Flux loss correction for isolated sources
        qover: Overlap correction for crowded peaks
        total: Combined correction (qcon * qover)
        qcon_sigma: Uncertainty on qcon (if computed)
        qover_sigma: Uncertainty on qover (if computed)
    """
    qcon: float
    qover: float
    total: float
    qcon_sigma: Optional[float] = None
    qover_sigma: Optional[float] = None

    @property
    def total_sigma(self) -> Optional[float]:
        """Combined uncertainty (simple quadrature sum)."""
        if self.qcon_sigma is None or self.qover_sigma is None:
            return None
        # Relative errors add in quadrature
        rel_qcon = self.qcon_sigma / self.qcon if self.qcon != 0 else 0
        rel_qover = self.qover_sigma / self.qover if self.qover != 0 else 0
        rel_total = np.sqrt(rel_qcon**2 + rel_qover**2)
        return rel_total * self.total


def compute_corrections(
    cut_ratio: float,
    dist_stat: float,
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
    dist_stat_sigma: float = 0.0,
) -> CorrectionFactors:
    """
    Compute all correction factors for Fourier filtering.

    Args:
        cut_ratio: Ratio of filter cut length to source size
        dist_stat: Distance statistic for overlap
        kernel: Filter kernel type
        dist_stat_sigma: Uncertainty on distance statistic

    Returns:
        CorrectionFactors with qcon, qover, and combined values
    """
    qcon = flux_loss_correction(cut_ratio, kernel)
    qover = overlap_correction(dist_stat, kernel)

    # Compute uncertainties if dist_stat_sigma provided
    qcon_sigma = None  # Not implemented for flux loss
    qover_sigma = None
    if dist_stat_sigma > 0:
        qover_sigma = overlap_correction_sigma(dist_stat, dist_stat_sigma, kernel)

    return CorrectionFactors(
        qcon=float(qcon) if np.isscalar(qcon) else qcon,
        qover=float(qover) if np.isscalar(qover) else qover,
        total=float(qcon * qover),
        qcon_sigma=qcon_sigma,
        qover_sigma=qover_sigma,
    )

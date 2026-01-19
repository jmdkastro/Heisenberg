"""
Diffuse Emission Fraction Calculation.

Calculates the fraction of diffuse (large-scale) vs compact (small-scale)
emission in astronomical images using Fourier filtering.

The diffuse fraction is defined as:
    diffuse_frac = 1 - compact_flux / total_flux

where compact flux is measured after highpass filtering.

References:
    IDL: fourier_diffuse_fraction.pro
    Kruijssen et al. (2018), MNRAS 479, 1866
"""

import numpy as np
from typing import Optional, Tuple, Union
from dataclasses import dataclass

from .filters import FilterKernel, FilterPass
from .filter_tool import apply_fourier_filter
from .corrections import flux_loss_correction, overlap_correction


@dataclass
class DiffuseFractionResult:
    """
    Result of diffuse fraction calculation.

    Attributes:
        flux_frac: Fraction of flux in compact component
        diffuse_frac: Fraction of flux in diffuse component (1 - flux_frac)
        flux_frac_err: Uncertainty on flux_frac (if computed)
        diffuse_frac_err: Uncertainty on diffuse_frac (if computed)
        total_flux: Total flux in original image
        compact_flux: Flux in highpass-filtered (compact) component
        diffuse_flux: Flux in lowpass-filtered (diffuse) component
    """
    flux_frac: float
    diffuse_frac: float
    flux_frac_err: Optional[float] = None
    diffuse_frac_err: Optional[float] = None
    total_flux: Optional[float] = None
    compact_flux: Optional[float] = None
    diffuse_flux: Optional[float] = None


def calculate_diffuse_fraction(
    image: np.ndarray,
    cut_length: float,
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
    order: int = 2,
    pass_type: Union[FilterPass, str] = FilterPass.HIGHPASS,
    mask: Optional[np.ndarray] = None,
    positive_only: bool = True,
) -> DiffuseFractionResult:
    """
    Calculate the diffuse emission fraction in an image.

    Args:
        image: 2D input image array
        cut_length: Filter cutoff length in pixels
        kernel: Filter kernel type ('gaussian', 'butterworth', 'ideal')
        order: Butterworth order (only used for Butterworth filters)
        pass_type: Filter pass type ('high' for compact, 'low' for diffuse)
        mask: Optional mask array (1=valid, 0=masked)
        positive_only: If True, only sum non-negative flux (default True)

    Returns:
        DiffuseFractionResult with flux fractions and components

    Notes:
        - With pass_type='high', filtered image contains compact emission
        - diffuse_frac = 1 - (compact_flux / total_flux)
        - NaN pixels are automatically excluded from flux sums
    """
    # Apply filter
    result = apply_fourier_filter(
        image,
        cut_length=cut_length,
        kernel=kernel,
        pass_type=pass_type,
        order=order,
        compute_residual=True,
    )

    filtered = result.filtered
    residual = result.residual

    # Apply mask if provided
    if mask is not None:
        image = np.where(mask > 0, image, np.nan)
        filtered = np.where(mask > 0, filtered, np.nan)
        if residual is not None:
            residual = np.where(mask > 0, residual, np.nan)

    # Calculate fluxes
    if positive_only:
        # Only sum non-negative values
        total_flux = np.nansum(np.where(image >= 0, image, 0))
        filtered_flux = np.nansum(np.where(filtered >= 0, filtered, 0))
    else:
        total_flux = np.nansum(image)
        filtered_flux = np.nansum(filtered)

    # Avoid division by zero
    if total_flux == 0:
        return DiffuseFractionResult(
            flux_frac=0.0,
            diffuse_frac=1.0,
            total_flux=0.0,
            compact_flux=0.0,
            diffuse_flux=0.0,
        )

    # Calculate flux fraction based on pass type
    if isinstance(pass_type, str):
        pass_type = FilterPass(pass_type.lower())

    if pass_type == FilterPass.HIGHPASS:
        # Filtered = compact emission
        compact_flux = filtered_flux
        diffuse_flux = total_flux - compact_flux
        flux_frac = compact_flux / total_flux
    else:
        # Filtered = diffuse emission
        diffuse_flux = filtered_flux
        compact_flux = total_flux - diffuse_flux
        flux_frac = compact_flux / total_flux

    diffuse_frac = 1.0 - flux_frac

    return DiffuseFractionResult(
        flux_frac=flux_frac,
        diffuse_frac=diffuse_frac,
        total_flux=total_flux,
        compact_flux=compact_flux,
        diffuse_flux=diffuse_flux,
    )


def calculate_diffuse_fraction_with_errors(
    image: np.ndarray,
    cut_length: float,
    cut_length_errmin: float,
    cut_length_errmax: float,
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
    order: int = 2,
    mask: Optional[np.ndarray] = None,
    positive_only: bool = True,
) -> DiffuseFractionResult:
    """
    Calculate diffuse fraction with error propagation from cut length uncertainty.

    Runs three filter calculations:
    - Nominal cut length
    - Cut length + errmax
    - Cut length - errmin

    Error bounds are derived from the difference in flux fractions.

    Args:
        image: 2D input image array
        cut_length: Nominal filter cutoff length in pixels
        cut_length_errmin: Lower error bound on cut length
        cut_length_errmax: Upper error bound on cut length
        kernel: Filter kernel type
        order: Butterworth order
        mask: Optional mask array
        positive_only: If True, only sum non-negative flux

    Returns:
        DiffuseFractionResult with flux fractions and uncertainties
    """
    # Nominal calculation
    result_mid = calculate_diffuse_fraction(
        image, cut_length, kernel, order,
        mask=mask, positive_only=positive_only
    )

    # Upper cut length (more diffuse)
    cut_max = cut_length + cut_length_errmax
    result_max = calculate_diffuse_fraction(
        image, cut_max, kernel, order,
        mask=mask, positive_only=positive_only
    )

    # Lower cut length (less diffuse)
    cut_min = max(cut_length - cut_length_errmin, 1.0)  # Avoid zero/negative
    result_min = calculate_diffuse_fraction(
        image, cut_min, kernel, order,
        mask=mask, positive_only=positive_only
    )

    # Error propagation
    # When cut_length increases, more flux passes lowpass -> diffuse_frac increases
    flux_frac_errmax = abs(result_max.flux_frac - result_mid.flux_frac)
    flux_frac_errmin = abs(result_mid.flux_frac - result_min.flux_frac)

    # Average error for simplicity
    flux_frac_err = (flux_frac_errmax + flux_frac_errmin) / 2.0
    diffuse_frac_err = flux_frac_err  # Same magnitude, opposite sign

    return DiffuseFractionResult(
        flux_frac=result_mid.flux_frac,
        diffuse_frac=result_mid.diffuse_frac,
        flux_frac_err=flux_frac_err,
        diffuse_frac_err=diffuse_frac_err,
        total_flux=result_mid.total_flux,
        compact_flux=result_mid.compact_flux,
        diffuse_flux=result_mid.diffuse_flux,
    )


def calculate_corrected_diffuse_fraction(
    image: np.ndarray,
    cut_length: float,
    cut_ratio: float,
    dist_stat: float,
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
    order: int = 2,
    mask: Optional[np.ndarray] = None,
) -> Tuple[DiffuseFractionResult, float, float]:
    """
    Calculate diffuse fraction with flux loss and overlap corrections.

    Args:
        image: 2D input image array
        cut_length: Filter cutoff length in pixels
        cut_ratio: Ratio of cut length to source size (for flux loss correction)
        dist_stat: Distance statistic (for overlap correction)
        kernel: Filter kernel type
        order: Butterworth order
        mask: Optional mask array

    Returns:
        Tuple of (result, qcon, qover) where:
            result: DiffuseFractionResult with corrected fractions
            qcon: Flux loss correction factor
            qover: Overlap correction factor

    Notes:
        Corrected compact flux = raw compact flux / (qcon * qover)
    """
    # Get raw result
    result = calculate_diffuse_fraction(
        image, cut_length, kernel, order, mask=mask
    )

    # Get correction factors
    qcon = float(flux_loss_correction(cut_ratio, kernel))
    qover = float(overlap_correction(dist_stat, kernel))
    total_correction = qcon * qover

    if total_correction == 0:
        return result, qcon, qover

    # Apply corrections to compact flux
    corrected_compact = result.compact_flux / total_correction

    # Recompute fractions with corrected compact flux
    # Note: total flux doesn't change, only attribution between compact/diffuse
    if result.total_flux > 0:
        corrected_flux_frac = corrected_compact / result.total_flux
        # Clamp to valid range
        corrected_flux_frac = max(0.0, min(1.0, corrected_flux_frac))
        corrected_diffuse_frac = 1.0 - corrected_flux_frac
    else:
        corrected_flux_frac = 0.0
        corrected_diffuse_frac = 1.0

    corrected_result = DiffuseFractionResult(
        flux_frac=corrected_flux_frac,
        diffuse_frac=corrected_diffuse_frac,
        total_flux=result.total_flux,
        compact_flux=corrected_compact,
        diffuse_flux=result.total_flux - corrected_compact,
    )

    return corrected_result, qcon, qover

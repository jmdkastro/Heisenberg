"""
Fourier filtering module for Heisenberg.

Contains Butterworth, Gaussian, and ideal filter implementations,
diffuse fraction calculation, and correction algorithms.

Main entry points:
    apply_fourier_filter: Apply frequency domain filter to image
    calculate_diffuse_fraction: Compute diffuse/compact emission split
    compute_corrections: Get flux loss and overlap corrections

Example:
    >>> from heisenberg.fourier import apply_fourier_filter, FilterKernel
    >>> result = apply_fourier_filter(image, cut_length=50, kernel='gaussian')
"""

from .filters import (
    FilterKernel,
    FilterPass,
    lowpass_gaussian,
    lowpass_butterworth,
    lowpass_ideal,
    highpass_gaussian,
    highpass_butterworth,
    highpass_ideal,
    get_filter,
    create_frequency_grid,
)

from .filter_tool import (
    FilterResult,
    apply_fourier_filter,
    compute_power_spectrum,
    compute_power_fraction,
    filter_to_resolution,
    extract_compact_emission,
)

from .corrections import (
    symmetric_sigmoidal,
    asymmetric_sigmoidal,
    flux_loss_correction,
    overlap_correction,
    overlap_correction_sigma,
    CorrectionFactors,
    compute_corrections,
    FLUX_LOSS_PARAMS_GAUSSIAN,
    OVERLAP_PARAMS_GAUSSIAN,
    OVERLAP_COVARIANCE_GAUSSIAN,
)

from .diffuse import (
    DiffuseFractionResult,
    calculate_diffuse_fraction,
    calculate_diffuse_fraction_with_errors,
    calculate_corrected_diffuse_fraction,
)

__all__ = [
    # Filter types
    "FilterKernel",
    "FilterPass",
    # Individual filters
    "lowpass_gaussian",
    "lowpass_butterworth",
    "lowpass_ideal",
    "highpass_gaussian",
    "highpass_butterworth",
    "highpass_ideal",
    "get_filter",
    "create_frequency_grid",
    # Filter tool
    "FilterResult",
    "apply_fourier_filter",
    "compute_power_spectrum",
    "compute_power_fraction",
    "filter_to_resolution",
    "extract_compact_emission",
    # Corrections
    "symmetric_sigmoidal",
    "asymmetric_sigmoidal",
    "flux_loss_correction",
    "overlap_correction",
    "overlap_correction_sigma",
    "CorrectionFactors",
    "compute_corrections",
    "FLUX_LOSS_PARAMS_GAUSSIAN",
    "OVERLAP_PARAMS_GAUSSIAN",
    "OVERLAP_COVARIANCE_GAUSSIAN",
    # Diffuse fraction
    "DiffuseFractionResult",
    "calculate_diffuse_fraction",
    "calculate_diffuse_fraction_with_errors",
    "calculate_corrected_diffuse_fraction",
]

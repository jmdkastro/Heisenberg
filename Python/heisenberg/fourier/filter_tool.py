"""
Fourier Filter Tool.

Core pipeline for applying Fourier domain filters to images.
Handles FFT, filter application, and inverse FFT with proper
NaN handling and precision control.

References:
    IDL: fourier_filter_tool.pro
"""

import numpy as np
from typing import Optional, Tuple, Union
from dataclasses import dataclass

from .filters import (
    FilterKernel,
    FilterPass,
    get_filter,
    create_frequency_grid,
)


@dataclass
class FilterResult:
    """
    Result of Fourier filtering operation.

    Attributes:
        filtered: The filtered image
        residual: Original minus filtered (None if not computed)
        taper: The filter taper values used
    """
    filtered: np.ndarray
    residual: Optional[np.ndarray] = None
    taper: Optional[np.ndarray] = None


def apply_fourier_filter(
    image: np.ndarray,
    cut_length: float,
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
    pass_type: Union[FilterPass, str] = FilterPass.HIGHPASS,
    order: int = 2,
    compute_residual: bool = False,
    return_taper: bool = False,
) -> FilterResult:
    """
    Apply a Fourier domain filter to an image.

    Args:
        image: 2D input image array
        cut_length: Cutoff length in pixels (converted to frequency as 1/cut_length)
        kernel: Filter kernel type ('gaussian', 'butterworth', 'ideal')
        pass_type: Filter pass type ('low' or 'high')
        order: Butterworth order (only used for Butterworth filters)
        compute_residual: If True, also compute original - filtered
        return_taper: If True, include the filter taper in result

    Returns:
        FilterResult with filtered image and optional residual/taper

    Notes:
        - NaN pixels are temporarily set to zero for FFT, then restored
        - Only the real part of the inverse FFT is returned
    """
    if image.ndim != 2:
        raise ValueError(f"Image must be 2D, got {image.ndim}D")

    # Track NaN locations
    nan_mask = np.isnan(image)
    has_nan = np.any(nan_mask)

    # Replace NaN with zero for FFT
    if has_nan:
        image_clean = np.where(nan_mask, 0.0, image)
    else:
        image_clean = image

    # Compute FFT
    fft_arr = np.fft.fft2(image_clean)

    # Create frequency grid and compute cutoff frequency
    freq_dist = create_frequency_grid(image.shape)
    cut_freq = 1.0 / cut_length

    # Get filter taper values
    taper = get_filter(kernel, pass_type, freq_dist, cut_freq, order)

    # Apply filter in frequency domain
    fft_filtered = fft_arr * taper

    # Inverse FFT and take real part
    filtered = np.real(np.fft.ifft2(fft_filtered))

    # Restore NaN pixels
    if has_nan:
        filtered = np.where(nan_mask, np.nan, filtered)

    # Compute residual if requested
    residual = None
    if compute_residual:
        residual = image - filtered

    return FilterResult(
        filtered=filtered,
        residual=residual,
        taper=taper if return_taper else None,
    )


def compute_power_spectrum(image: np.ndarray) -> np.ndarray:
    """
    Compute the power spectrum of an image.

    Args:
        image: 2D input image array

    Returns:
        Power spectrum (|FFT|²)
    """
    # Handle NaN
    image_clean = np.where(np.isnan(image), 0.0, image)
    fft_arr = np.fft.fft2(image_clean)
    return np.abs(fft_arr)**2


def compute_power_fraction(
    image: np.ndarray,
    cut_lengths: np.ndarray,
) -> np.ndarray:
    """
    Compute fraction of power retained after highpass filtering.

    For each cut length, computes what fraction of total power
    remains after applying an ideal highpass filter.

    Args:
        image: 2D input image array
        cut_lengths: Array of cut lengths in pixels

    Returns:
        Array of power fractions (same length as cut_lengths)

    Notes:
        Uses ideal highpass filter for clean power separation.
    """
    # Compute power spectrum
    power = compute_power_spectrum(image)
    total_power = np.sum(power)

    if total_power == 0:
        return np.zeros(len(cut_lengths))

    # Create frequency grid
    freq_dist = create_frequency_grid(image.shape)

    # Compute power fraction for each cut length
    fractions = np.zeros(len(cut_lengths))

    for i, cut_length in enumerate(cut_lengths):
        cut_freq = 1.0 / cut_length
        # Use ideal highpass for clean separation
        taper = get_filter('ideal', 'high', freq_dist, cut_freq)
        partial_power = np.sum(power * taper)
        fractions[i] = partial_power / total_power

    return fractions


def filter_to_resolution(
    image: np.ndarray,
    target_length: float,
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
    order: int = 2,
) -> np.ndarray:
    """
    Lowpass filter image to remove structure smaller than target length.

    Convenience function for smoothing in Fourier space.

    Args:
        image: 2D input image array
        target_length: Target resolution in pixels (structures smaller
            than this will be attenuated)
        kernel: Filter kernel type
        order: Butterworth order (if applicable)

    Returns:
        Filtered image with small-scale structure removed
    """
    result = apply_fourier_filter(
        image,
        cut_length=target_length,
        kernel=kernel,
        pass_type=FilterPass.LOWPASS,
        order=order,
    )
    return result.filtered


def extract_compact_emission(
    image: np.ndarray,
    cut_length: float,
    kernel: Union[FilterKernel, str] = FilterKernel.GAUSSIAN,
    order: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Separate compact and diffuse emission using highpass filtering.

    Args:
        image: 2D input image array
        cut_length: Separation scale in pixels
        kernel: Filter kernel type
        order: Butterworth order (if applicable)

    Returns:
        Tuple of (compact_emission, diffuse_emission)

    Notes:
        compact = highpass filtered (small-scale structure)
        diffuse = original - compact (large-scale structure)
    """
    result = apply_fourier_filter(
        image,
        cut_length=cut_length,
        kernel=kernel,
        pass_type=FilterPass.HIGHPASS,
        order=order,
        compute_residual=True,
    )
    return result.filtered, result.residual

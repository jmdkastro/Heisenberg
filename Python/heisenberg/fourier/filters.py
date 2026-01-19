"""
Fourier Domain Filters.

Implements lowpass and highpass filters for spatial frequency filtering:
- Gaussian: Smooth exponential rolloff
- Butterworth: Parameterized rolloff sharpness
- Ideal: Sharp binary cutoff

These filters are applied in Fourier space to separate diffuse and
compact emission components in astronomical images.

References:
    IDL: fourier_lowpass_*.pro, fourier_highpass_*.pro
"""

import numpy as np
from enum import Enum
from typing import Union


class FilterKernel(Enum):
    """Available filter kernel types."""
    GAUSSIAN = "gaussian"
    BUTTERWORTH = "butterworth"
    IDEAL = "ideal"


class FilterPass(Enum):
    """Filter pass type."""
    LOWPASS = "low"
    HIGHPASS = "high"


def lowpass_gaussian(
    freq_dist: np.ndarray,
    cut_freq: float,
) -> np.ndarray:
    """
    Gaussian lowpass filter.

    Smooth exponential rolloff that gradually attenuates high frequencies.

    Args:
        freq_dist: Array of frequency distances from FFT origin
        cut_freq: Cutoff frequency (1/cut_length)

    Returns:
        Filter taper values in range [0, 1]

    Formula:
        taper = exp(-(freq_dist² / (2 × cut_freq²)))
    """
    return np.exp(-(freq_dist**2) / (2.0 * cut_freq**2))


def lowpass_butterworth(
    freq_dist: np.ndarray,
    cut_freq: float,
    order: int = 2,
) -> np.ndarray:
    """
    Butterworth lowpass filter.

    Parameterized rolloff where higher orders give sharper cutoffs.

    Args:
        freq_dist: Array of frequency distances from FFT origin
        cut_freq: Cutoff frequency (1/cut_length)
        order: Butterworth order (positive integer, default 2)

    Returns:
        Filter taper values in range [0, 1]

    Formula:
        taper = 1 / (1 + (freq_dist/cut_freq)^(2*order))

    Raises:
        ValueError: If order is not a positive integer
    """
    if not isinstance(order, int) or order < 1:
        raise ValueError(f"Butterworth order must be positive integer, got {order}")

    # Avoid division by zero at zero frequency
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = freq_dist / cut_freq
        taper = 1.0 / (1.0 + ratio**(2 * order))
        # Handle zero frequency (ratio=0 -> taper=1)
        taper = np.where(np.isfinite(taper), taper, 1.0)

    return taper


def lowpass_ideal(
    freq_dist: np.ndarray,
    cut_freq: float,
    equal: bool = False,
) -> np.ndarray:
    """
    Ideal lowpass filter (sharp binary cutoff).

    Args:
        freq_dist: Array of frequency distances from FFT origin
        cut_freq: Cutoff frequency (1/cut_length)
        equal: If True, use <= instead of < for cutoff

    Returns:
        Filter taper values (0 or 1)

    Note:
        May cause ringing artifacts (Gibbs phenomenon) due to sharp cutoff.
    """
    if equal:
        return np.where(freq_dist <= cut_freq, 1.0, 0.0)
    else:
        return np.where(freq_dist < cut_freq, 1.0, 0.0)


def highpass_gaussian(
    freq_dist: np.ndarray,
    cut_freq: float,
) -> np.ndarray:
    """
    Gaussian highpass filter.

    Complement of Gaussian lowpass: 1 - lowpass_gaussian.

    Args:
        freq_dist: Array of frequency distances from FFT origin
        cut_freq: Cutoff frequency (1/cut_length)

    Returns:
        Filter taper values in range [0, 1]

    Formula:
        taper = 1 - exp(-(freq_dist² / (2 × cut_freq²)))
    """
    return 1.0 - lowpass_gaussian(freq_dist, cut_freq)


def highpass_butterworth(
    freq_dist: np.ndarray,
    cut_freq: float,
    order: int = 2,
) -> np.ndarray:
    """
    Butterworth highpass filter.

    Note: This is NOT simply 1 - lowpass. The formula uses inverted ratio.

    Args:
        freq_dist: Array of frequency distances from FFT origin
        cut_freq: Cutoff frequency (1/cut_length)
        order: Butterworth order (positive integer, default 2)

    Returns:
        Filter taper values in range [0, 1]

    Formula:
        taper = 1 / (1 + (cut_freq/freq_dist)^(2*order))

    Raises:
        ValueError: If order is not a positive integer
    """
    if not isinstance(order, int) or order < 1:
        raise ValueError(f"Butterworth order must be positive integer, got {order}")

    # Avoid division by zero at zero frequency
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = cut_freq / freq_dist
        taper = 1.0 / (1.0 + ratio**(2 * order))
        # Handle zero frequency (ratio=inf -> taper=0)
        taper = np.where(np.isfinite(taper), taper, 0.0)

    return taper


def highpass_ideal(
    freq_dist: np.ndarray,
    cut_freq: float,
    equal: bool = False,
) -> np.ndarray:
    """
    Ideal highpass filter (sharp binary cutoff).

    Args:
        freq_dist: Array of frequency distances from FFT origin
        cut_freq: Cutoff frequency (1/cut_length)
        equal: If True, use >= instead of > for cutoff

    Returns:
        Filter taper values (0 or 1)
    """
    if equal:
        return np.where(freq_dist >= cut_freq, 1.0, 0.0)
    else:
        return np.where(freq_dist > cut_freq, 1.0, 0.0)


def get_filter(
    kernel: Union[FilterKernel, str],
    pass_type: Union[FilterPass, str],
    freq_dist: np.ndarray,
    cut_freq: float,
    order: int = 2,
) -> np.ndarray:
    """
    Get filter taper values for specified kernel and pass type.

    Convenience function that dispatches to the appropriate filter function.

    Args:
        kernel: Filter kernel type ('gaussian', 'butterworth', 'ideal')
        pass_type: Filter pass type ('low' or 'high')
        freq_dist: Array of frequency distances from FFT origin
        cut_freq: Cutoff frequency (1/cut_length)
        order: Butterworth order (only used for Butterworth filters)

    Returns:
        Filter taper values

    Raises:
        ValueError: If kernel or pass_type is invalid
    """
    # Convert strings to enums
    if isinstance(kernel, str):
        kernel = FilterKernel(kernel.lower())
    if isinstance(pass_type, str):
        pass_type = FilterPass(pass_type.lower())

    if pass_type == FilterPass.LOWPASS:
        if kernel == FilterKernel.GAUSSIAN:
            return lowpass_gaussian(freq_dist, cut_freq)
        elif kernel == FilterKernel.BUTTERWORTH:
            return lowpass_butterworth(freq_dist, cut_freq, order)
        elif kernel == FilterKernel.IDEAL:
            return lowpass_ideal(freq_dist, cut_freq)
    else:  # HIGHPASS
        if kernel == FilterKernel.GAUSSIAN:
            return highpass_gaussian(freq_dist, cut_freq)
        elif kernel == FilterKernel.BUTTERWORTH:
            return highpass_butterworth(freq_dist, cut_freq, order)
        elif kernel == FilterKernel.IDEAL:
            return highpass_ideal(freq_dist, cut_freq)

    raise ValueError(f"Invalid kernel: {kernel}")


def create_frequency_grid(shape: tuple) -> np.ndarray:
    """
    Create a 2D frequency distance grid for FFT filtering.

    The grid contains the distance of each frequency from the origin
    (DC component), normalized such that the Nyquist frequency is 0.5.

    Args:
        shape: Shape of the image (ny, nx)

    Returns:
        2D array of frequency distances, same shape as input

    Notes:
        - Uses numpy's fftfreq which handles even/odd dimensions correctly
        - Origin (DC) is at [0, 0] after fftshift
        - Frequencies are in cycles per pixel
    """
    ny, nx = shape

    # Create 1D frequency vectors
    freq_x = np.fft.fftfreq(nx)
    freq_y = np.fft.fftfreq(ny)

    # Create 2D meshgrid
    fx, fy = np.meshgrid(freq_x, freq_y)

    # Compute distance from origin
    freq_dist = np.sqrt(fx**2 + fy**2)

    return freq_dist

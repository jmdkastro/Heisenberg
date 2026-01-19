"""
Image convolution utilities for Heisenberg.

This module provides functions for convolving astronomical images with
Gaussian or tophat kernels, handling NaN values appropriately.
"""

import numpy as np
from typing import Tuple, Optional, Union
from scipy.ndimage import convolve
from scipy.signal import fftconvolve

from astropy.convolution import Gaussian2DKernel, convolve_fft


def psf_gaussian(
    fwhm: float,
    size: Optional[int] = None,
    normalize: bool = True
) -> np.ndarray:
    """
    Generate a 2D Gaussian PSF kernel.

    Args:
        fwhm: Full width at half maximum in pixels
        size: Kernel size in pixels (default: 5 * FWHM, odd)
        normalize: If True, normalize kernel to sum to 1

    Returns:
        2D Gaussian kernel array
    """
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    if size is None:
        size = int(5 * fwhm)
        # Ensure odd size
        size = 2 * (size // 2) + 1
        size = max(size, 3)

    # Use astropy's Gaussian2DKernel for accurate normalization
    kernel = Gaussian2DKernel(x_stddev=sigma, y_stddev=sigma, x_size=size, y_size=size)

    if normalize:
        return kernel.array / np.sum(kernel.array)
    return kernel.array


def psf_tophat(width: float) -> np.ndarray:
    """
    Generate a 2D tophat (circular pillbox) PSF kernel.

    This accurately computes the fractional coverage of each pixel by the
    circular tophat, accounting for partial pixel coverage at the edges.

    Args:
        width: Diameter of the tophat in pixels

    Returns:
        2D tophat kernel array, normalized to sum to 1
    """
    npix = 2 * int(round(width / 2.0)) + 1  # Always odd
    frac = np.zeros((npix, npix))
    psf = np.zeros((npix, npix))

    middle = 0.5 * (npix - 1.0)
    radius = 0.5 * width

    for i in range(npix):
        for j in range(npix):
            # Pixel corner coordinates relative to center
            x1 = i - 0.5 - middle
            x2 = i + 0.5 - middle
            y1 = j - 0.5 - middle
            y2 = j + 0.5 - middle

            # Corner coordinates and radii
            x = np.array([x1, x1, x2, x2])
            y = np.array([y1, y2, y1, y2])
            rad = np.sqrt(x**2 + y**2)

            # Sort by radius
            sortidx = np.argsort(rad)
            x = x[sortidx]
            y = y[sortidx]
            rad = rad[sortidx]

            # Count corners inside circle
            nin = np.sum(rad <= radius)

            if nin == 0:
                # Entirely outside - check for edge cases
                frac[i, j] = 0.0
                if max(abs(x)) / max(abs(y)) > 1:
                    rmin = min(abs(x))
                else:
                    rmin = min(abs(y))

                # Central row/column edge case
                if (i == int(middle) or j == int(middle)) and rmin <= radius:
                    frac[i, j] = (radius**2 * np.arccos(rmin / radius) -
                                  rmin * np.sqrt(2 * radius * (radius - rmin) -
                                               (radius - rmin)**2))

                # Central pixel case
                if i == int(middle) and j == int(middle):
                    if radius >= 0.5:
                        frac[i, j] = (np.pi * radius**2 -
                                     4 * (radius**2 * np.arccos(0.5 / radius) -
                                          0.5 * np.sqrt(2 * radius * (radius - 0.5) -
                                                       (radius - 0.5)**2)))
                    else:
                        frac[i, j] = np.pi * radius**2

            elif nin == 4:
                # Entirely inside
                frac[i, j] = 1.0

            else:
                # Partial coverage - compute geometric intersection
                if nin == 1:
                    lx = np.sqrt(radius**2 - y[0]**2) - abs(x[0])
                    ly = np.sqrt(radius**2 - x[0]**2) - abs(y[0])
                    a = np.sqrt(lx**2 + ly**2)
                    area = 0.5 * lx * ly

                elif nin == 2:
                    if x[0] == x[1]:
                        ccst, cvar = x, y
                    else:
                        ccst, cvar = y, x
                    l1 = np.sqrt(radius**2 - cvar[0]**2) - abs(ccst[0])
                    l2 = np.sqrt(radius**2 - cvar[1]**2) - abs(ccst[1])
                    dl = l1 - l2
                    a = np.sqrt(1 + dl**2)
                    rmax = max(abs(np.concatenate([x, y])))
                    if (i == int(middle) or j == int(middle)) and rmax <= radius:
                        minarea = (radius**2 * np.arccos(rmax / radius) -
                                  rmax * np.sqrt(2 * radius * (radius - rmax) -
                                               (radius - rmax)**2))
                    else:
                        minarea = 0.0
                    area = l2 + 0.5 * dl - minarea

                elif nin == 3:
                    lx = abs(x[3]) - np.sqrt(radius**2 - y[3]**2)
                    ly = abs(y[3]) - np.sqrt(radius**2 - x[3]**2)
                    a = np.sqrt(lx**2 + ly**2)
                    area = 1.0 - 0.5 * lx * ly

                # Add circular segment area
                rcs = 0.5 * np.sqrt(4 * radius**2 - a**2)
                hcs = radius - rcs
                areacs = radius**2 * np.arccos(rcs / radius) - rcs * np.sqrt(2 * radius * hcs - hcs**2)
                frac[i, j] = area + areacs

            # Normalize by circle area
            psf[i, j] = frac[i, j] / (np.pi * radius**2)

    return psf


def convolve_image(
    image: np.ndarray,
    kernel: np.ndarray,
    preserve_nan: bool = True,
    use_fft: bool = True
) -> np.ndarray:
    """
    Convolve an image with a kernel, properly handling NaN values.

    Args:
        image: Input 2D image
        kernel: Convolution kernel
        preserve_nan: If True, restore NaN positions after convolution
        use_fft: If True, use FFT convolution (faster for large kernels)

    Returns:
        Convolved image
    """
    # Handle NaN values
    nan_mask = np.isnan(image)
    has_nans = np.any(nan_mask)

    if has_nans:
        # Replace NaNs with zeros for convolution
        image_clean = np.where(nan_mask, 0.0, image)
    else:
        image_clean = image

    # Perform convolution
    if use_fft:
        result = convolve_fft(image_clean, kernel, normalize_kernel=True,
                             nan_treatment='fill', fill_value=0.0,
                             preserve_nan=False)
    else:
        result = convolve(image_clean, kernel, mode='constant', cval=0.0)
        result = result / np.sum(kernel)

    # Restore NaN positions
    if preserve_nan and has_nans:
        result[nan_mask] = np.nan

    return result


def smooth_to_resolution(
    image: np.ndarray,
    current_fwhm: float,
    target_fwhm: float,
    pixel_scale: float,
    tophat: bool = False,
    preserve_nan: bool = True
) -> np.ndarray:
    """
    Smooth an image to a target resolution.

    Args:
        image: Input 2D image
        current_fwhm: Current resolution FWHM in arcseconds
        target_fwhm: Target resolution FWHM in arcseconds
        pixel_scale: Pixel scale in arcseconds per pixel
        tophat: If True, use tophat kernel; if False, use Gaussian
        preserve_nan: If True, restore NaN positions after convolution

    Returns:
        Smoothed image at target resolution

    Raises:
        ValueError: If target resolution is finer than current resolution
    """
    if target_fwhm < current_fwhm * 0.95 and not tophat:
        raise ValueError(
            f"Target resolution ({target_fwhm:.1f}\") is finer than "
            f"current resolution ({current_fwhm:.1f}\")"
        )

    if tophat:
        # Tophat kernel with diameter = target FWHM
        kernel_width = target_fwhm / pixel_scale  # pixels
        kernel = psf_tophat(kernel_width)
    else:
        # Gaussian kernel to convolve from current to target resolution
        convolve_fwhm = np.sqrt(target_fwhm**2 - current_fwhm**2)
        kernel_fwhm = convolve_fwhm / pixel_scale  # pixels
        kernel = psf_gaussian(kernel_fwhm)

    return convolve_image(image, kernel, preserve_nan=preserve_nan)


def regrid_image(
    image: np.ndarray,
    scale_factor: float,
    order: int = 3
) -> np.ndarray:
    """
    Regrid (resample) an image to a coarser pixel scale.

    Args:
        image: Input 2D image
        scale_factor: Factor by which to increase pixel size (>1 = coarser)
        order: Interpolation order (0=nearest, 1=linear, 3=cubic)

    Returns:
        Regridded image with fewer pixels
    """
    from scipy.ndimage import zoom

    if scale_factor <= 1.0:
        return image.copy()

    # Handle NaN values
    nan_mask = np.isnan(image)
    has_nans = np.any(nan_mask)

    if has_nans:
        image_clean = np.where(nan_mask, 0.0, image)
        nan_weight = (~nan_mask).astype(float)
    else:
        image_clean = image

    # Zoom to new size
    zoom_factor = 1.0 / scale_factor
    result = zoom(image_clean, zoom_factor, order=order)

    # Handle NaN propagation
    if has_nans:
        nan_result = zoom(nan_weight, zoom_factor, order=order)
        result = np.where(nan_result < 0.5, np.nan, result)

    return result

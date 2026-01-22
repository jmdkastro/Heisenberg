"""
High-Level Peak Detection.

Provides a user-friendly interface for peak detection using the clumpfind
algorithm. Handles contour level generation, sensitivity filtering, and
dual-map (star/gas) detection.

Based on the IDL peaks2d.pro and peak_find.pro implementations.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .clumpfind import clumpfind2d, ClumpfindResult
from .statistics import (
    PeakStatistics,
    compute_clump_statistics,
    compute_all_statistics,
)


@dataclass
class PeakDetectionConfig:
    """
    Configuration for peak detection.

    Attributes:
        npixmin: Minimum pixels per clump (default: 20)
        nsigma: Sensitivity threshold multiplier (default: 5.0)
        loglevels: Use logarithmic level spacing (default: True)
        logrange: Range in dex for log spacing (default: 2.0)
        logspacing: Interval between log levels in dex (default: 0.5)
        nlevels: Number of contour levels (calculated externally as
            logrange / logspacing + 1 for log, or nlinlevel for linear)
        flux_weighted: Use flux-weighted positions (default: False)
    """
    npixmin: int = 20
    nsigma: float = 5.0
    loglevels: bool = True
    logrange: float = 2.0
    logspacing: float = 0.5
    nlevels: int = 5  # Default: logrange/logspacing + 1 = 2.0/0.5 + 1 = 5
    flux_weighted: bool = False


@dataclass
class DetectedPeak:
    """
    A detected peak with position, flux, and optional full statistics.

    Attributes:
        x: X position (column)
        y: Y position (row)
        total_flux: Total integrated flux
        npix: Number of pixels
        peak_flux: Peak flux value
        stats: Full PeakStatistics if requested
    """
    x: float
    y: float
    total_flux: float
    npix: int
    peak_flux: float = 0.0
    stats: Optional[PeakStatistics] = None


def generate_levels(
    image: np.ndarray,
    nlevels: int,
    loglevels: bool = True,
    logrange: float = 2.0,
    logspacing: float = 0.5,
) -> np.ndarray:
    """
    Generate contour levels for clumpfind.

    This function is equivalent to the level generation in IDL peak_find.pro.
    It receives nlevels as a parameter (calculated externally, like in IDL
    tuningfork.pro as: nlevels = logrange / logspacing + 1 for log levels,
    or nlevels = nlinlevel for linear levels).

    Args:
        image: Input image to determine level range from
        nlevels: Number of contour levels (pre-calculated)
        loglevels: If True, use logarithmic spacing; if False, linear
        logrange: For log spacing, the range in dex below the maximum
        logspacing: For log spacing, the interval between levels in dex

    Returns:
        Array of contour levels from lowest to highest

    Notes:
        For logarithmic spacing (from IDL peak_find.pro):
            maxval = max(alog10(peakid_x), /nan)
            maxlevel = (floor(maxval/logspacing_x) - 1) * logspacing_x
            levels = 10^(maxlevel - logrange_x + dindgen(nlevels_x)/(nlevels_x-1) * logrange_x)

        For linear spacing (from IDL peak_find.pro):
            maxval = max(peakid_x, /nan)
            minval = min(peakid_x, /nan)
            levels = minval + (maxval - minval) * dindgen(nlevels_x) / (nlevels_x - 1)

    References:
        IDL: peak_find.pro
    """
    # Get valid (non-NaN) values
    valid = image[~np.isnan(image)]
    if len(valid) == 0:
        raise ValueError("Image contains no valid (non-NaN) values")

    maxval = np.max(valid)

    if maxval <= 0:
        raise ValueError("Image maximum must be positive for level generation")

    if loglevels:
        # Logarithmic spacing from IDL peak_find.pro
        # IDL: maxval = max(alog10(peakid_x), /nan)
        log_max = np.log10(maxval)

        # IDL: maxlevel = (floor(maxval/logspacing_x) - 1) * logspacing_x
        # (Note: IDL maxval is already log10 of image)
        maxlevel = (np.floor(log_max / logspacing) - 1) * logspacing

        # IDL: levels = 10.^(maxlevel - logrange_x + dindgen(nlevels_x)/(nlevels_x-1) * logrange_x)
        log_levels = maxlevel - logrange + np.arange(nlevels) / (nlevels - 1) * logrange
        levels = 10**log_levels
    else:
        # Linear spacing from IDL peak_find.pro
        # IDL: maxval = max(peakid_x, /nan)
        # IDL: minval = min(peakid_x, /nan)
        minval = np.min(valid)

        # IDL: levels = minval + (maxval - minval) * dindgen(nlevels_x) / (nlevels_x - 1)
        levels = minval + (maxval - minval) * np.arange(nlevels) / (nlevels - 1)

    return levels


def generate_contour_levels(
    image: np.ndarray,
    config: PeakDetectionConfig,
) -> np.ndarray:
    """
    Generate contour levels for clumpfind based on configuration.

    Wrapper around generate_levels that extracts parameters from config.

    Args:
        image: Input image to determine level range
        config: PeakDetectionConfig with level parameters

    Returns:
        Array of contour levels from lowest to highest
    """
    return generate_levels(
        image,
        nlevels=config.nlevels,
        loglevels=config.loglevels,
        logrange=config.logrange,
        logspacing=config.logspacing,
    )


def find_peaks(
    image: np.ndarray,
    config: Optional[PeakDetectionConfig] = None,
    sensitivity: Optional[np.ndarray] = None,
    offset: float = 0.0,
    include_stats: bool = True,
) -> List[DetectedPeak]:
    """
    Find peaks in an image using the clumpfind algorithm.

    This is the main entry point for peak detection. It:
    1. Generates contour levels based on configuration
    2. Runs the clumpfind algorithm
    3. Computes statistics for each clump
    4. Filters by sensitivity threshold
    5. Returns peaks sorted by total flux (descending)

    Args:
        image: 2D input image array. NaN values are masked.
        config: PeakDetectionConfig with detection parameters.
            If None, uses default configuration.
        sensitivity: Optional 2D sensitivity map (same shape as image).
            Used for filtering: peak_flux > nsigma * sensitivity + offset
        offset: Offset for sensitivity filtering (default: 0.0)
        include_stats: If True, include full PeakStatistics in result

    Returns:
        List of DetectedPeak objects sorted by total flux (descending)

    Example:
        >>> config = PeakDetectionConfig(npixmin=20, nsigma=5.0)
        >>> peaks = find_peaks(image, config)
        >>> print(f"Found {len(peaks)} peaks")
    """
    if config is None:
        config = PeakDetectionConfig()

    # Generate contour levels
    levels = generate_contour_levels(image, config)

    # Run clumpfind
    result = clumpfind2d(image, levels, npixmin=config.npixmin)

    if not result.clumps:
        return []

    # Compute statistics for each clump
    all_stats = compute_all_statistics(image, result, config.flux_weighted)

    # Build detected peaks
    peaks = []
    for clump, stats in zip(result.clumps, all_stats):
        peak = DetectedPeak(
            x=stats.x,
            y=stats.y,
            total_flux=stats.total_flux,
            npix=stats.npix,
            peak_flux=stats.peak_flux,
            stats=stats if include_stats else None,
        )
        peaks.append(peak)

    # Apply sensitivity filtering if provided
    if sensitivity is not None:
        peaks = _filter_by_sensitivity(
            peaks, image, sensitivity, config.nsigma, offset
        )

    # Sort by total flux (descending)
    peaks.sort(key=lambda p: p.total_flux, reverse=True)

    return peaks


def _filter_by_sensitivity(
    peaks: List[DetectedPeak],
    image: np.ndarray,
    sensitivity: np.ndarray,
    nsigma: float,
    offset: float
) -> List[DetectedPeak]:
    """
    Filter peaks by sensitivity threshold.

    A peak is kept if: peak_flux > nsigma * sensitivity_at_peak + offset

    Args:
        peaks: List of DetectedPeak objects
        image: Original image (for shape)
        sensitivity: Sensitivity map
        nsigma: Multiplier for sensitivity threshold
        offset: Offset added to threshold

    Returns:
        Filtered list of peaks
    """
    if sensitivity.shape != image.shape:
        raise ValueError(
            f"Sensitivity shape {sensitivity.shape} must match "
            f"image shape {image.shape}"
        )

    filtered = []
    for peak in peaks:
        # Get sensitivity at peak position
        row = int(round(peak.y))
        col = int(round(peak.x))

        # Bounds check
        row = max(0, min(row, image.shape[0] - 1))
        col = max(0, min(col, image.shape[1] - 1))

        sens_at_peak = sensitivity[row, col]
        threshold = nsigma * sens_at_peak + offset

        if peak.peak_flux > threshold:
            filtered.append(peak)

    return filtered


def find_peaks_dual(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    star_config: Optional[PeakDetectionConfig] = None,
    gas_config: Optional[PeakDetectionConfig] = None,
    star_sensitivity: Optional[np.ndarray] = None,
    gas_sensitivity: Optional[np.ndarray] = None,
    star_offset: float = 0.0,
    gas_offset: float = 0.0,
) -> Tuple[List[DetectedPeak], List[DetectedPeak]]:
    """
    Find peaks in both star and gas maps.

    Convenience function for the typical Heisenberg workflow where
    peaks need to be identified in both stellar and gas tracer maps.

    Args:
        star_image: Stellar tracer image
        gas_image: Gas tracer image
        star_config: Configuration for star peak detection
        gas_config: Configuration for gas peak detection
        star_sensitivity: Sensitivity map for star image
        gas_sensitivity: Sensitivity map for gas image
        star_offset: Offset for star sensitivity filtering
        gas_offset: Offset for gas sensitivity filtering

    Returns:
        Tuple of (star_peaks, gas_peaks)

    Example:
        >>> star_peaks, gas_peaks = find_peaks_dual(star_map, gas_map)
        >>> print(f"Star: {len(star_peaks)}, Gas: {len(gas_peaks)}")
    """
    star_peaks = find_peaks(
        star_image,
        config=star_config,
        sensitivity=star_sensitivity,
        offset=star_offset,
    )

    gas_peaks = find_peaks(
        gas_image,
        config=gas_config,
        sensitivity=gas_sensitivity,
        offset=gas_offset,
    )

    return star_peaks, gas_peaks


def peaks_to_array(peaks: List[DetectedPeak]) -> np.ndarray:
    """
    Convert list of DetectedPeak to a numpy array.

    Returns array with columns: [x, y, total_flux, npix]

    Args:
        peaks: List of DetectedPeak objects

    Returns:
        2D array with shape (n_peaks, 4)
    """
    if not peaks:
        return np.empty((0, 4))

    return np.array([
        [p.x, p.y, p.total_flux, p.npix]
        for p in peaks
    ])


def peaks_to_coords(peaks: List[DetectedPeak]) -> np.ndarray:
    """
    Extract peak coordinates as (x, y) array.

    Args:
        peaks: List of DetectedPeak objects

    Returns:
        2D array with shape (n_peaks, 2) containing [x, y] coordinates
    """
    if not peaks:
        return np.empty((0, 2))

    return np.array([[p.x, p.y] for p in peaks])

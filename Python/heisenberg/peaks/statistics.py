"""
Peak Statistics Calculation.

Computes detailed statistics for identified clumps, including position,
flux, size metrics (FWHM), and boundary detection.

Based on the IDL clstats2d.pro implementation.

References:
    Williams, de Geus, & Blitz 1994, ApJ, 428, 693
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Optional

from .clumpfind import Clump, ClumpfindResult


# Conversion factor from sigma to FWHM for Gaussian
SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))  # ~2.355


@dataclass
class PeakStatistics:
    """
    Statistics for a single identified peak/clump.

    Attributes:
        id: Clump identifier
        x: X position (column) - peak or flux-weighted centroid
        y: Y position (row) - peak or flux-weighted centroid
        peak_flux: Maximum flux value in the clump
        total_flux: Total integrated flux (sum of all pixels)
        fwhm_x: Full width at half maximum in x direction
        fwhm_y: Full width at half maximum in y direction
        radius: Equivalent circular radius = sqrt(npix / pi)
        npix: Number of pixels in the clump
        on_edge: True if clump touches image boundary
    """
    id: int
    x: float
    y: float
    peak_flux: float
    total_flux: float
    fwhm_x: float
    fwhm_y: float
    radius: float
    npix: int
    on_edge: bool = False


def compute_clump_statistics(
    image: np.ndarray,
    clump: Clump,
    flux_weighted: bool = False
) -> PeakStatistics:
    """
    Compute statistics for a single clump.

    Args:
        image: The original 2D image array
        clump: The Clump object to compute statistics for
        flux_weighted: If True, compute flux-weighted centroid position.
            If False, use peak pixel position.

    Returns:
        PeakStatistics object with computed values

    Notes:
        FWHM is computed from the flux-weighted second moment:
            variance_i = sum(i^2 * flux) / sum(flux) - ibar^2
            FWHM = 2.355 * sqrt(variance)
    """
    # Extract pixel coordinates and values
    rows = np.array([p[0] for p in clump.pixels])
    cols = np.array([p[1] for p in clump.pixels])
    fluxes = image[rows, cols]

    # Handle any NaN values in the flux
    valid = ~np.isnan(fluxes)
    if not np.any(valid):
        # All NaN - return minimal statistics
        return PeakStatistics(
            id=clump.id,
            x=float(clump.peak_position[1]),
            y=float(clump.peak_position[0]),
            peak_flux=np.nan,
            total_flux=np.nan,
            fwhm_x=np.nan,
            fwhm_y=np.nan,
            radius=np.sqrt(clump.npix / np.pi),
            npix=clump.npix,
            on_edge=_check_on_edge(rows, cols, image.shape)
        )

    rows = rows[valid]
    cols = cols[valid]
    fluxes = fluxes[valid]

    # Total and peak flux
    total_flux = float(np.sum(fluxes))
    peak_flux = float(np.max(fluxes))

    # Position
    if flux_weighted and total_flux > 0:
        # Flux-weighted centroid
        x = float(np.sum(cols * fluxes) / total_flux)
        y = float(np.sum(rows * fluxes) / total_flux)
    else:
        # Peak position
        x = float(clump.peak_position[1])
        y = float(clump.peak_position[0])

    # Compute FWHM from variance
    fwhm_x, fwhm_y = _compute_fwhm(rows, cols, fluxes, y, x)

    # Equivalent radius
    radius = np.sqrt(clump.npix / np.pi)

    # Check if on edge
    on_edge = _check_on_edge(rows, cols, image.shape)

    return PeakStatistics(
        id=clump.id,
        x=x,
        y=y,
        peak_flux=peak_flux,
        total_flux=total_flux,
        fwhm_x=fwhm_x,
        fwhm_y=fwhm_y,
        radius=radius,
        npix=clump.npix,
        on_edge=on_edge
    )


def _compute_fwhm(
    rows: np.ndarray,
    cols: np.ndarray,
    fluxes: np.ndarray,
    center_y: float,
    center_x: float
) -> tuple:
    """
    Compute FWHM in x and y directions from flux-weighted variance.

    The FWHM is computed as 2.355 * sqrt(variance), where variance
    is the flux-weighted second moment about the centroid.

    Args:
        rows: Row indices of pixels
        cols: Column indices of pixels
        fluxes: Flux values at each pixel
        center_y: Y position of centroid
        center_x: X position of centroid

    Returns:
        Tuple of (fwhm_x, fwhm_y)
    """
    total_flux = np.sum(fluxes)

    if total_flux <= 0:
        return (0.0, 0.0)

    # Flux-weighted second moment (variance)
    # var_x = sum((x - x_center)^2 * flux) / sum(flux)
    var_x = np.sum((cols - center_x)**2 * fluxes) / total_flux
    var_y = np.sum((rows - center_y)**2 * fluxes) / total_flux

    # FWHM = 2.355 * sigma, where sigma = sqrt(variance)
    fwhm_x = SIGMA_TO_FWHM * np.sqrt(max(var_x, 0))
    fwhm_y = SIGMA_TO_FWHM * np.sqrt(max(var_y, 0))

    return (float(fwhm_x), float(fwhm_y))


def _check_on_edge(
    rows: np.ndarray,
    cols: np.ndarray,
    image_shape: tuple
) -> bool:
    """
    Check if any pixel in the clump is on the image boundary.

    Args:
        rows: Row indices of clump pixels
        cols: Column indices of clump pixels
        image_shape: Shape of the image (nrows, ncols)

    Returns:
        True if any pixel touches the edge
    """
    nrows, ncols = image_shape

    # Check boundaries
    if np.any(rows == 0) or np.any(rows == nrows - 1):
        return True
    if np.any(cols == 0) or np.any(cols == ncols - 1):
        return True

    return False


def compute_all_statistics(
    image: np.ndarray,
    result: ClumpfindResult,
    flux_weighted: bool = False
) -> List[PeakStatistics]:
    """
    Compute statistics for all clumps in a ClumpfindResult.

    Args:
        image: The original 2D image array
        result: ClumpfindResult from clumpfind2d
        flux_weighted: If True, compute flux-weighted centroid positions

    Returns:
        List of PeakStatistics objects, one per clump
    """
    return [
        compute_clump_statistics(image, clump, flux_weighted)
        for clump in result.clumps
    ]


def statistics_to_array(stats_list: List[PeakStatistics]) -> np.ndarray:
    """
    Convert a list of PeakStatistics to a numpy array.

    Returns array with columns:
    [id, x, y, peak_flux, total_flux, fwhm_x, fwhm_y, radius, npix]

    Args:
        stats_list: List of PeakStatistics objects

    Returns:
        2D numpy array with shape (n_peaks, 9)
    """
    if not stats_list:
        return np.empty((0, 9))

    data = np.array([
        [
            s.id,
            s.x,
            s.y,
            s.peak_flux,
            s.total_flux,
            s.fwhm_x,
            s.fwhm_y,
            s.radius,
            s.npix
        ]
        for s in stats_list
    ])

    return data

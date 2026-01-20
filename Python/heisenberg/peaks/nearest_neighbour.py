"""
Nearest Neighbour Statistics for Peak Analysis.

Calculates median relative nearest-neighbour distances for peaks,
used for statistical significance analysis.

References:
    IDL: med_peak_relative_nearest_neighbour_dist_calc.pro
"""

from typing import Optional, Tuple, List
from dataclasses import dataclass

import numpy as np
from scipy.spatial.distance import cdist


@dataclass
class NearestNeighbourResult:
    """
    Result from nearest-neighbour distance calculation.

    Attributes:
        dist_med: Median nearest-neighbour distance normalized by FWHM
        dist_med_sigma: Error on dist_med
        dist_val: Median nearest-neighbour distance in pixels
        dist_sigma: Error on dist_val
        fwhm_val: Mean FWHM of peaks
        fwhm_sigma: Error on fwhm_val
        n_peaks: Number of peaks used
    """
    dist_med: float
    dist_med_sigma: float
    dist_val: float
    dist_sigma: float
    fwhm_val: float
    fwhm_sigma: float
    n_peaks: int


def med_peak_relative_nearest_neighbour_dist(
    x: np.ndarray,
    y: np.ndarray,
    fwhm_x: np.ndarray,
    fwhm_y: np.ndarray,
    mask: Optional[np.ndarray] = None,
    n_bootstrap: int = 10000,
    seed: Optional[int] = None,
) -> NearestNeighbourResult:
    """
    Calculate median relative nearest-neighbour distance for peaks.

    This is a Python translation of IDL med_peak_relative_nearest_neighbour_dist_calc.pro.

    The relative nearest-neighbour distance is the median of the minimum
    distances between peaks, divided by the mean FWHM of the peaks.
    Bootstrap resampling is used to estimate uncertainties.

    Args:
        x: X coordinates of peaks (pixels)
        y: Y coordinates of peaks (pixels)
        fwhm_x: X-direction FWHM of each peak (pixels)
        fwhm_y: Y-direction FWHM of each peak (pixels)
        mask: Optional mask array; if provided, only include peaks where
              mask[y, x] == 1
        n_bootstrap: Number of bootstrap iterations for error estimation
        seed: Random seed for reproducibility

    Returns:
        NearestNeighbourResult with distances and uncertainties
    """
    rng = np.random.default_rng(seed)

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    fwhm_x = np.asarray(fwhm_x, dtype=float)
    fwhm_y = np.asarray(fwhm_y, dtype=float)

    # Apply mask if provided
    if mask is not None:
        # Find peaks inside mask
        x_int = np.round(x).astype(int)
        y_int = np.round(y).astype(int)

        # Clip to valid indices
        x_int = np.clip(x_int, 0, mask.shape[1] - 1)
        y_int = np.clip(y_int, 0, mask.shape[0] - 1)

        valid = mask[y_int, x_int] == 1
        x = x[valid]
        y = y[valid]
        fwhm_x = fwhm_x[valid]
        fwhm_y = fwhm_y[valid]

    n_peaks = len(x)

    if n_peaks < 2:
        return NearestNeighbourResult(
            dist_med=np.nan,
            dist_med_sigma=np.nan,
            dist_val=np.nan,
            dist_sigma=np.nan,
            fwhm_val=np.nan,
            fwhm_sigma=np.nan,
            n_peaks=n_peaks,
        )

    # Calculate mean FWHM (average of x and y for each peak, then mean over peaks)
    fwhm_vec = (fwhm_x + fwhm_y) / 2
    fwhm_all = np.concatenate([fwhm_x, fwhm_y])
    fwhm_val = np.mean(fwhm_all)

    # Bootstrap FWHM to get error
    fboot_quant = np.zeros(n_bootstrap)
    for b in range(n_bootstrap):
        indices = rng.integers(0, len(fwhm_all), size=len(fwhm_all))
        fboot_quant[b] = np.mean(fwhm_all[indices])

    bfwhm_val = np.median(fboot_quant)
    bfwhm_lq = np.percentile(fboot_quant, 16)
    bfwhm_uq = np.percentile(fboot_quant, 84)
    fwhm_sigma = np.mean(np.abs([bfwhm_lq - bfwhm_val, bfwhm_uq - bfwhm_val]))

    # Calculate pairwise distances
    coords = np.column_stack([x, y])
    dist_matrix = cdist(coords, coords, metric='euclidean')

    # Set diagonal to NaN (distance to self)
    np.fill_diagonal(dist_matrix, np.nan)

    # Get minimum distance for each peak (nearest neighbour)
    min_dist_vec = np.nanmin(dist_matrix, axis=1)

    # Median of nearest-neighbour distances
    dist_val = np.median(min_dist_vec)

    # Bootstrap min_dist_vec to get error
    n_dists = len(min_dist_vec)
    boot_quant = np.zeros(n_bootstrap)
    for b in range(n_bootstrap):
        indices = rng.integers(0, n_dists, size=n_dists)
        boot_quant[b] = np.median(min_dist_vec[indices])

    dist_lq = np.percentile(boot_quant, 16)
    dist_uq = np.percentile(boot_quant, 84)
    dist_sigma = np.mean(np.abs([dist_lq - dist_val, dist_uq - dist_val]))

    # Calculate final normalized value and error
    dist_med = dist_val / fwhm_val
    dist_med_sigma = dist_med * np.sqrt(
        (dist_sigma / dist_val) ** 2 + (fwhm_sigma / fwhm_val) ** 2
    )

    return NearestNeighbourResult(
        dist_med=dist_med,
        dist_med_sigma=dist_med_sigma,
        dist_val=dist_val,
        dist_sigma=dist_sigma,
        fwhm_val=fwhm_val,
        fwhm_sigma=fwhm_sigma,
        n_peaks=n_peaks,
    )


def nearest_neighbour_from_peaks(
    peaks: List,
    mask: Optional[np.ndarray] = None,
    n_bootstrap: int = 10000,
    seed: Optional[int] = None,
) -> NearestNeighbourResult:
    """
    Calculate nearest-neighbour statistics from a list of DetectedPeak objects.

    Args:
        peaks: List of DetectedPeak objects
        mask: Optional mask array
        n_bootstrap: Number of bootstrap iterations
        seed: Random seed

    Returns:
        NearestNeighbourResult with distances and uncertainties
    """
    if not peaks:
        return NearestNeighbourResult(
            dist_med=np.nan,
            dist_med_sigma=np.nan,
            dist_val=np.nan,
            dist_sigma=np.nan,
            fwhm_val=np.nan,
            fwhm_sigma=np.nan,
            n_peaks=0,
        )

    x = np.array([p.x for p in peaks])
    y = np.array([p.y for p in peaks])

    # Get FWHM from stats if available
    fwhm_x = np.array([
        p.stats.fwhm_x if hasattr(p, 'stats') and p.stats is not None else 1.0
        for p in peaks
    ])
    fwhm_y = np.array([
        p.stats.fwhm_y if hasattr(p, 'stats') and p.stats is not None else 1.0
        for p in peaks
    ])

    return med_peak_relative_nearest_neighbour_dist(
        x, y, fwhm_x, fwhm_y,
        mask=mask,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )

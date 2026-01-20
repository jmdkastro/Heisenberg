"""
Peak identification module for Heisenberg.

Contains the CLFIND2D algorithm for 2D clump finding, peak detection
wrappers, and peak statistics calculation.

Main entry points:
    find_peaks: High-level peak detection function
    find_peaks_dual: Detect peaks in both star and gas maps
    clumpfind2d: Low-level clumpfind algorithm

Example:
    >>> from heisenberg.peaks import find_peaks, PeakDetectionConfig
    >>> config = PeakDetectionConfig(npixmin=20, nsigma=5.0)
    >>> peaks = find_peaks(image, config)
"""

from .clumpfind import (
    Clump,
    ClumpfindResult,
    clumpfind2d,
    generate_levels,
)

from .statistics import (
    PeakStatistics,
    compute_clump_statistics,
    compute_all_statistics,
    statistics_to_array,
    SIGMA_TO_FWHM,
)

from .detection import (
    PeakDetectionConfig,
    DetectedPeak,
    find_peaks,
    find_peaks_dual,
    generate_contour_levels,
    peaks_to_array,
    peaks_to_coords,
)

from .nearest_neighbour import (
    NearestNeighbourResult,
    med_peak_relative_nearest_neighbour_dist,
    nearest_neighbour_from_peaks,
)

__all__ = [
    # Clumpfind
    "Clump",
    "ClumpfindResult",
    "clumpfind2d",
    "generate_levels",
    # Statistics
    "PeakStatistics",
    "compute_clump_statistics",
    "compute_all_statistics",
    "statistics_to_array",
    "SIGMA_TO_FWHM",
    # Detection
    "PeakDetectionConfig",
    "DetectedPeak",
    "find_peaks",
    "find_peaks_dual",
    "generate_contour_levels",
    "peaks_to_array",
    "peaks_to_coords",
    # Nearest neighbour statistics
    "NearestNeighbourResult",
    "med_peak_relative_nearest_neighbour_dist",
    "nearest_neighbour_from_peaks",
]

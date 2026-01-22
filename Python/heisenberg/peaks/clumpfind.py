"""
CLFIND2D - Hierarchical Clump Finding Algorithm.

Implementation of the Williams, de Geus, & Blitz (1994, ApJ, 428, 693) algorithm
for identifying clumps in 2D images using a hierarchical contour-based approach.

The algorithm works by:
1. Processing contour levels from highest to lowest
2. At each level, identifying connected regions (8-connectivity)
3. Creating new clumps when peaks are found that aren't already assigned
4. When clumps merge at lower levels, assigning pixels to the nearest peak

References:
    Williams, de Geus, & Blitz 1994, ApJ, 428, 693
    IDL implementation: clfind2d.pro
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional
from scipy import ndimage


@dataclass
class Clump:
    """
    A single identified clump in the image.

    Attributes:
        id: Unique identifier for this clump (1-indexed)
        peak_position: (row, col) position of the peak pixel
        peak_value: Flux value at the peak pixel
        pixels: List of (row, col) tuples for all pixels in the clump
        npix: Number of pixels in the clump
    """
    id: int
    peak_position: Tuple[int, int]
    peak_value: float
    pixels: List[Tuple[int, int]] = field(default_factory=list)
    npix: int = 0

    def __post_init__(self):
        if self.npix == 0 and self.pixels:
            self.npix = len(self.pixels)


@dataclass
class ClumpfindResult:
    """
    Result of the clumpfind algorithm.

    Attributes:
        clumps: List of identified Clump objects, sorted by decreasing peak flux
        assignment_map: 2D array with same shape as input, where each pixel
            contains the clump ID it belongs to (0 = unassigned)
        levels: The contour levels used for clump identification
    """
    clumps: List[Clump]
    assignment_map: np.ndarray
    levels: np.ndarray


def clumpfind2d(
    image: np.ndarray,
    levels: np.ndarray,
    npixmin: int = 20,
) -> ClumpfindResult:
    """
    Hierarchical clump identification algorithm.

    Identifies clumps in a 2D image by processing contour levels from highest
    to lowest. At each level, connected regions are found and either assigned
    to existing clumps or used to create new ones.

    Args:
        image: 2D input image array. NaN values are treated as masked.
        levels: Array of contour levels, will be processed from highest to lowest.
        npixmin: Minimum number of pixels required for a valid clump.
            Clumps with fewer pixels are rejected. Default is 20.

    Returns:
        ClumpfindResult containing the identified clumps, assignment map,
        and the levels used.

    Notes:
        - Uses 8-connectivity (including diagonals) for region identification
        - When clumps merge at lower contour levels, pixels are assigned to
          the clump with the nearest peak (Euclidean distance)
        - Clumps are sorted by decreasing peak flux in the final result
    """
    if image.ndim != 2:
        raise ValueError(f"Image must be 2D, got {image.ndim}D")

    # Sort levels from highest to lowest
    levels = np.sort(levels)[::-1]

    # Initialize assignment map (0 = unassigned)
    assignment_map = np.zeros(image.shape, dtype=np.int32)

    # Track clump peaks: {clump_id: (row, col)}
    clump_peaks: Dict[int, Tuple[int, int]] = {}
    clump_peak_values: Dict[int, float] = {}

    # Next available clump ID
    next_clump_id = 1

    # 8-connectivity structure for scipy.ndimage.label
    connectivity = np.ones((3, 3), dtype=np.int32)

    # Create mask for valid pixels (not NaN)
    valid_mask = ~np.isnan(image)

    # Process levels from highest to lowest
    for level in levels:
        # Find pixels above this level
        above_level = (image >= level) & valid_mask

        # Find connected regions at this level
        labeled, n_regions = ndimage.label(above_level, structure=connectivity)

        # Process each region
        for region_id in range(1, n_regions + 1):
            region_mask = labeled == region_id
            region_pixels = np.where(region_mask)

            if len(region_pixels[0]) == 0:
                continue

            # Find which clumps this region overlaps with
            overlapping_clumps = set(assignment_map[region_mask]) - {0}

            if len(overlapping_clumps) == 0:
                # No overlap - this is a new clump
                # Find the peak pixel in this region
                region_values = image[region_mask]
                peak_idx = np.argmax(region_values)
                peak_row = region_pixels[0][peak_idx]
                peak_col = region_pixels[1][peak_idx]
                peak_value = region_values[peak_idx]

                # Create new clump
                clump_id = next_clump_id
                next_clump_id += 1

                clump_peaks[clump_id] = (peak_row, peak_col)
                clump_peak_values[clump_id] = peak_value

                # Assign all pixels in this region to the new clump
                assignment_map[region_mask] = clump_id

            elif len(overlapping_clumps) == 1:
                # Single clump - extend it
                clump_id = overlapping_clumps.pop()
                # Assign unassigned pixels in this region to the clump
                unassigned = region_mask & (assignment_map == 0)
                assignment_map[unassigned] = clump_id

            else:
                # Multiple clumps merge - assign pixels to nearest peak
                _assign_merged_pixels(
                    region_mask,
                    assignment_map,
                    clump_peaks,
                    overlapping_clumps
                )

    # Filter clumps by minimum pixel count and build final result
    clumps = _filter_and_build_clumps(
        assignment_map, clump_peaks, clump_peak_values, npixmin
    )

    return ClumpfindResult(
        clumps=clumps,
        assignment_map=assignment_map,
        levels=levels
    )


def _assign_merged_pixels(
    region_mask: np.ndarray,
    assignment_map: np.ndarray,
    clump_peaks: Dict[int, Tuple[int, int]],
    overlapping_clumps: set
) -> None:
    """
    Assign pixels in a merged region to the nearest clump peak.

    When multiple clumps merge at a lower contour level, each pixel
    is assigned to the clump whose peak is closest (Euclidean distance).

    Args:
        region_mask: Boolean mask of the merged region
        assignment_map: Assignment map to update in place
        clump_peaks: Dictionary mapping clump IDs to peak positions
        overlapping_clumps: Set of clump IDs that overlap this region
    """
    # Get all pixels that need assignment (currently unassigned in this region)
    unassigned = region_mask & (assignment_map == 0)
    unassigned_coords = np.where(unassigned)

    if len(unassigned_coords[0]) == 0:
        return

    # For each unassigned pixel, find the nearest clump peak
    for i in range(len(unassigned_coords[0])):
        row = unassigned_coords[0][i]
        col = unassigned_coords[1][i]

        min_dist_sq = np.inf
        nearest_clump = None

        for clump_id in overlapping_clumps:
            peak_row, peak_col = clump_peaks[clump_id]
            dist_sq = (row - peak_row)**2 + (col - peak_col)**2

            if dist_sq < min_dist_sq:
                min_dist_sq = dist_sq
                nearest_clump = clump_id

        if nearest_clump is not None:
            assignment_map[row, col] = nearest_clump


def _filter_and_build_clumps(
    assignment_map: np.ndarray,
    clump_peaks: Dict[int, Tuple[int, int]],
    clump_peak_values: Dict[int, float],
    npixmin: int
) -> List[Clump]:
    """
    Filter clumps by minimum pixel count and build final Clump objects.

    Args:
        assignment_map: The pixel assignment map
        clump_peaks: Dictionary mapping clump IDs to peak positions
        clump_peak_values: Dictionary mapping clump IDs to peak values
        npixmin: Minimum number of pixels required

    Returns:
        List of Clump objects sorted by decreasing peak flux
    """
    clumps = []

    for clump_id in clump_peaks:
        pixels = np.where(assignment_map == clump_id)
        npix = len(pixels[0])

        if npix < npixmin:
            # Mark as rejected (set to 0)
            assignment_map[assignment_map == clump_id] = 0
            continue

        pixel_list = list(zip(pixels[0].tolist(), pixels[1].tolist()))

        clump = Clump(
            id=clump_id,
            peak_position=clump_peaks[clump_id],
            peak_value=clump_peak_values[clump_id],
            pixels=pixel_list,
            npix=npix
        )
        clumps.append(clump)

    # Sort by decreasing peak flux
    clumps.sort(key=lambda c: c.peak_value, reverse=True)

    # Renumber clumps sequentially and update assignment map
    new_assignment = np.zeros_like(assignment_map)
    for new_id, clump in enumerate(clumps, start=1):
        old_id = clump.id
        new_assignment[assignment_map == old_id] = new_id
        clump.id = new_id

    # Update the assignment map in place
    assignment_map[:] = new_assignment

    return clumps



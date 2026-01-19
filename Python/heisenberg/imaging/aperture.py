"""
Aperture Flux Measurement.

Measures flux at peak positions across multiple aperture sizes using
pre-smoothed image cubes. This is the core photometry for the tuning
fork diagram.

The approach uses smoothed maps where each pixel already represents
the aperture-convolved flux, so flux at a peak position gives the
integrated flux within that aperture size.

References:
    IDL: tuningfork.pro lines 656-676
"""

import numpy as np
from typing import List, Optional, Tuple
from dataclasses import dataclass

from .convolution import smooth_to_resolution


@dataclass
class ApertureFluxResult:
    """
    Result of aperture flux measurement.

    Attributes:
        star_flux: Array of shape (n_apertures, n_peaks) with stellar flux
        gas_flux: Array of shape (n_apertures, n_peaks) with gas flux
        aperture_area_frac: Fractional area available at each peak/aperture
        apertures: Array of aperture sizes used (in pixels or physical units)
        peak_positions: Array of (x, y) peak positions
    """
    star_flux: np.ndarray
    gas_flux: np.ndarray
    aperture_area_frac: np.ndarray
    apertures: np.ndarray
    peak_positions: np.ndarray


def generate_aperture_sizes(
    lap_min: float,
    lap_max: float,
    n_apertures: int,
) -> np.ndarray:
    """
    Generate geometric series of aperture sizes.

    Args:
        lap_min: Minimum aperture size
        lap_max: Maximum aperture size
        n_apertures: Number of apertures

    Returns:
        Array of aperture sizes from lap_min to lap_max

    Notes:
        Uses geometric spacing: apertures = lap_min * (lap_max/lap_min)^(i/(n-1))
    """
    if n_apertures == 1:
        return np.array([lap_min])

    ratio = lap_max / lap_min
    indices = np.arange(n_apertures) / (n_apertures - 1)
    return lap_min * (ratio ** indices)


def create_smoothed_cube(
    image: np.ndarray,
    apertures: np.ndarray,
    pixel_scale: float,
    beam_fwhm: float = 0.0,
    use_tophat: bool = True,
) -> np.ndarray:
    """
    Create a cube of smoothed images at different aperture sizes.

    Args:
        image: 2D input image
        apertures: Array of aperture sizes (in same units as pixel_scale)
        pixel_scale: Pixel scale (e.g., pc/pixel)
        beam_fwhm: Current beam FWHM (0 if no existing beam)
        use_tophat: If True, use tophat kernel; if False, use Gaussian

    Returns:
        3D array of shape (n_apertures, ny, nx) with smoothed images

    Notes:
        For each aperture, the image is smoothed to that resolution.
        The smoothed value at a pixel represents the mean flux within
        an aperture of that size centered on that pixel.
    """
    n_apertures = len(apertures)
    ny, nx = image.shape
    cube = np.zeros((n_apertures, ny, nx))

    for i, aperture in enumerate(apertures):
        # Convert aperture to FWHM (aperture is diameter)
        target_fwhm = aperture

        if target_fwhm <= beam_fwhm:
            # Already at or below target resolution
            cube[i] = image.copy()
        else:
            cube[i] = smooth_to_resolution(
                image,
                current_fwhm=beam_fwhm,
                target_fwhm=target_fwhm,
                pixel_scale=pixel_scale,
                tophat=use_tophat,
            )

    return cube


def measure_aperture_flux(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    star_peaks: np.ndarray,
    gas_peaks: np.ndarray,
    apertures: np.ndarray,
    pixel_scale: float,
    star_beam_fwhm: float = 0.0,
    gas_beam_fwhm: float = 0.0,
    mask: Optional[np.ndarray] = None,
    use_tophat: bool = True,
) -> Tuple[ApertureFluxResult, ApertureFluxResult]:
    """
    Measure flux at peaks across multiple aperture sizes.

    Args:
        star_image: Stellar tracer map
        gas_image: Gas tracer map
        star_peaks: Array of stellar peak positions, shape (n_star, 2) with (x, y)
        gas_peaks: Array of gas peak positions, shape (n_gas, 2) with (x, y)
        apertures: Array of aperture sizes
        pixel_scale: Pixel scale (physical units per pixel)
        star_beam_fwhm: Current beam FWHM for stellar image
        gas_beam_fwhm: Current beam FWHM for gas image
        mask: Optional mask array (1=valid, 0=masked)
        use_tophat: Use tophat kernel for smoothing

    Returns:
        Tuple of (star_result, gas_result) ApertureFluxResult objects
        - star_result: Flux at stellar peaks
        - gas_result: Flux at gas peaks
    """
    n_apertures = len(apertures)
    n_star_peaks = len(star_peaks)
    n_gas_peaks = len(gas_peaks)

    # Create smoothed cubes
    smooth_star = create_smoothed_cube(
        star_image, apertures, pixel_scale, star_beam_fwhm, use_tophat
    )
    smooth_gas = create_smoothed_cube(
        gas_image, apertures, pixel_scale, gas_beam_fwhm, use_tophat
    )

    # Create smoothed mask cube if mask provided
    if mask is not None:
        smooth_mask = create_smoothed_cube(
            mask.astype(float), apertures, pixel_scale, 0.0, use_tophat
        )
    else:
        smooth_mask = np.ones((n_apertures,) + star_image.shape)

    # Measure flux at stellar peaks
    star_star_flux = np.zeros((n_apertures, n_star_peaks))
    star_gas_flux = np.zeros((n_apertures, n_star_peaks))
    star_area_frac = np.zeros((n_apertures, n_star_peaks))

    for j, (px, py) in enumerate(star_peaks):
        # Convert to integer pixel indices
        ix, iy = int(round(px)), int(round(py))
        # Bounds check
        if 0 <= iy < star_image.shape[0] and 0 <= ix < star_image.shape[1]:
            for i in range(n_apertures):
                star_star_flux[i, j] = smooth_star[i, iy, ix]
                star_gas_flux[i, j] = smooth_gas[i, iy, ix]
                star_area_frac[i, j] = smooth_mask[i, iy, ix]

    star_result = ApertureFluxResult(
        star_flux=star_star_flux,
        gas_flux=star_gas_flux,
        aperture_area_frac=star_area_frac,
        apertures=apertures.copy(),
        peak_positions=star_peaks.copy(),
    )

    # Measure flux at gas peaks
    gas_star_flux = np.zeros((n_apertures, n_gas_peaks))
    gas_gas_flux = np.zeros((n_apertures, n_gas_peaks))
    gas_area_frac = np.zeros((n_apertures, n_gas_peaks))

    for j, (px, py) in enumerate(gas_peaks):
        ix, iy = int(round(px)), int(round(py))
        if 0 <= iy < gas_image.shape[0] and 0 <= ix < gas_image.shape[1]:
            for i in range(n_apertures):
                gas_star_flux[i, j] = smooth_star[i, iy, ix]
                gas_gas_flux[i, j] = smooth_gas[i, iy, ix]
                gas_area_frac[i, j] = smooth_mask[i, iy, ix]

    gas_result = ApertureFluxResult(
        star_flux=gas_star_flux,
        gas_flux=gas_gas_flux,
        aperture_area_frac=gas_area_frac,
        apertures=apertures.copy(),
        peak_positions=gas_peaks.copy(),
    )

    return star_result, gas_result


def compute_flux_ratios(
    flux_result: ApertureFluxResult,
    min_area_frac: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute gas-to-stellar flux ratios from aperture measurements.

    Args:
        flux_result: ApertureFluxResult from measure_aperture_flux
        min_area_frac: Minimum fractional area to include peak

    Returns:
        Tuple of (flux_ratios, valid_mask) where:
            flux_ratios: Shape (n_apertures,) mean flux ratio at each aperture
            valid_mask: Shape (n_apertures, n_peaks) boolean mask of valid measurements
    """
    n_apertures, n_peaks = flux_result.star_flux.shape

    # Create validity mask based on area fraction and non-zero flux
    valid = (
        (flux_result.aperture_area_frac >= min_area_frac) &
        (flux_result.star_flux > 0) &
        (flux_result.gas_flux > 0) &
        np.isfinite(flux_result.star_flux) &
        np.isfinite(flux_result.gas_flux)
    )

    # Compute flux ratios
    flux_ratios = np.zeros(n_apertures)

    for i in range(n_apertures):
        if np.any(valid[i]):
            total_gas = np.sum(flux_result.gas_flux[i, valid[i]])
            total_star = np.sum(flux_result.star_flux[i, valid[i]])
            if total_star > 0:
                flux_ratios[i] = total_gas / total_star

    return flux_ratios, valid


def select_non_overlapping_peaks(
    peak_positions: np.ndarray,
    min_separation: float,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """
    Select a subset of peaks that don't overlap within given separation.

    Uses random sampling to select non-overlapping peaks.

    Args:
        peak_positions: Array of (x, y) positions, shape (n_peaks, 2)
        min_separation: Minimum distance between selected peaks
        rng: Random number generator (uses default if None)

    Returns:
        Boolean mask indicating which peaks are selected
    """
    if rng is None:
        rng = np.random.default_rng()

    n_peaks = len(peak_positions)
    if n_peaks == 0:
        return np.array([], dtype=bool)

    # Shuffle order for random selection
    order = rng.permutation(n_peaks)

    selected = np.zeros(n_peaks, dtype=bool)
    selected_positions = []

    for idx in order:
        pos = peak_positions[idx]

        # Check distance to all already selected peaks
        if len(selected_positions) == 0:
            selected[idx] = True
            selected_positions.append(pos)
        else:
            distances = np.sqrt(np.sum((np.array(selected_positions) - pos)**2, axis=1))
            if np.all(distances >= min_separation):
                selected[idx] = True
                selected_positions.append(pos)

    return selected


def monte_carlo_flux_ratios(
    flux_result: ApertureFluxResult,
    n_mc: int = 100,
    min_area_frac: float = 0.5,
    seed: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute flux ratios with Monte Carlo sampling of non-overlapping peaks.

    Args:
        flux_result: ApertureFluxResult from measure_aperture_flux
        n_mc: Number of Monte Carlo realizations
        min_area_frac: Minimum fractional area to include peak
        seed: Random seed for reproducibility

    Returns:
        Tuple of (mean_ratios, std_ratios, n_peaks_used) where:
            mean_ratios: Shape (n_apertures,) mean flux ratio
            std_ratios: Shape (n_apertures,) standard deviation
            n_peaks_used: Shape (n_apertures,) mean number of peaks per MC sample
    """
    rng = np.random.default_rng(seed)
    n_apertures = len(flux_result.apertures)

    all_ratios = np.zeros((n_apertures, n_mc))
    all_n_peaks = np.zeros((n_apertures, n_mc))

    for i_ap, aperture in enumerate(flux_result.apertures):
        for i_mc in range(n_mc):
            # Select non-overlapping peaks for this aperture size
            selected = select_non_overlapping_peaks(
                flux_result.peak_positions,
                min_separation=aperture,
                rng=rng,
            )

            # Apply validity mask
            valid = (
                selected &
                (flux_result.aperture_area_frac[i_ap] >= min_area_frac) &
                (flux_result.star_flux[i_ap] > 0) &
                np.isfinite(flux_result.star_flux[i_ap]) &
                np.isfinite(flux_result.gas_flux[i_ap])
            )

            if np.any(valid):
                total_gas = np.sum(flux_result.gas_flux[i_ap, valid])
                total_star = np.sum(flux_result.star_flux[i_ap, valid])
                if total_star > 0:
                    all_ratios[i_ap, i_mc] = total_gas / total_star
                all_n_peaks[i_ap, i_mc] = np.sum(valid)

    mean_ratios = np.mean(all_ratios, axis=1)
    std_ratios = np.std(all_ratios, axis=1)
    mean_n_peaks = np.mean(all_n_peaks, axis=1)

    return mean_ratios, std_ratios, mean_n_peaks

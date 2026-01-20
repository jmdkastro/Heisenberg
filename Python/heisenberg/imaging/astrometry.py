"""
Astrometry utilities for Heisenberg.

Provides functions for WCS comparison and astrometric calculations.
This is a Python translation of IDL astrometry_equal.pro and related functions.

References:
    IDL: astrometry_equal.pro, get_platescale.pro
"""

from pathlib import Path
from typing import Optional, Tuple, Union, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from astropy.wcs import WCS


def get_equinox(header: dict) -> Optional[float]:
    """
    Get equinox from FITS header.

    Args:
        header: FITS header dictionary

    Returns:
        Equinox value (e.g., 2000.0 for J2000) or None if not found
    """
    # Check various keywords
    for key in ['EQUINOX', 'EPOCH', 'EQUINX']:
        if key in header:
            return float(header[key])

    # Check RADESYS
    radesys = header.get('RADESYS', header.get('RADECSYS', ''))
    if radesys.upper() in ('FK5', 'ICRS', 'J2000'):
        return 2000.0
    elif radesys.upper() in ('FK4', 'B1950'):
        return 1950.0

    return None


def get_platescale(
    wcs: 'WCS',
) -> Tuple[float, float]:
    """
    Get plate scale from WCS in degrees per pixel.

    This is a Python translation of IDL get_platescale.pro.

    Args:
        wcs: WCS object

    Returns:
        Tuple of (cdelt_x, cdelt_y) in degrees per pixel
    """
    from astropy.wcs import WCS

    if wcs is None:
        raise ValueError("No WCS provided")

    # Try to get CDELT values
    try:
        if hasattr(wcs.wcs, 'cd') and wcs.wcs.cd is not None:
            # CD matrix - extract scale
            cd = wcs.wcs.cd
            cdelt_x = np.sqrt(cd[0, 0]**2 + cd[1, 0]**2)
            cdelt_y = np.sqrt(cd[0, 1]**2 + cd[1, 1]**2)
            # Check sign from diagonal
            if cd[0, 0] < 0:
                cdelt_x = -cdelt_x
            if cd[1, 1] < 0:
                cdelt_y = -cdelt_y
        elif hasattr(wcs.wcs, 'cdelt') and wcs.wcs.cdelt is not None:
            cdelt_x = wcs.wcs.cdelt[0]
            cdelt_y = wcs.wcs.cdelt[1]
        else:
            # Try pixel scale utility
            scale = wcs.proj_plane_pixel_scales()
            cdelt_x = -scale[0]  # Typically RA decreases with x
            cdelt_y = scale[1]
    except Exception:
        raise ValueError("Cannot determine plate scale from WCS")

    return cdelt_x, cdelt_y


def get_rotation(
    wcs: 'WCS',
) -> float:
    """
    Get rotation angle from WCS in degrees.

    Args:
        wcs: WCS object

    Returns:
        Rotation angle in degrees
    """
    if wcs is None:
        raise ValueError("No WCS provided")

    try:
        if hasattr(wcs.wcs, 'cd') and wcs.wcs.cd is not None:
            cd = wcs.wcs.cd
            # Rotation from CD matrix
            rotation = np.rad2deg(np.arctan2(cd[0, 1], cd[0, 0]))
        elif hasattr(wcs.wcs, 'crota') and wcs.wcs.crota is not None:
            rotation = wcs.wcs.crota[1] if len(wcs.wcs.crota) > 1 else wcs.wcs.crota[0]
        else:
            rotation = 0.0
    except Exception:
        rotation = 0.0

    return rotation


def astrometry_equal(
    image_a: np.ndarray,
    header_a: Union[dict, 'WCS'],
    image_b: np.ndarray,
    header_b: Union[dict, 'WCS'],
    tolerance: float = 5e-4,
    strict_equinox: bool = False,
) -> bool:
    """
    Check if two images have the same astrometry.

    This is a Python translation of IDL astrometry_equal.pro.

    Args:
        image_a: First image array
        header_a: FITS header or WCS object for first image
        image_b: Second image array
        header_b: FITS header or WCS object for second image
        tolerance: Allowable tolerance in degrees between pixel positions
                  (default: 5e-4 ~ 1.8 arcsec)
        strict_equinox: If True, require exact equinox match

    Returns:
        True if astrometry is equal within tolerance, False otherwise

    Example:
        >>> from heisenberg.io.fits import read_fits
        >>> img_a, _, hdr_a = read_fits('image1.fits')
        >>> img_b, _, hdr_b = read_fits('image2.fits')
        >>> if astrometry_equal(img_a, hdr_a, img_b, hdr_b):
        ...     print("Images have matching astrometry")
    """
    from astropy.wcs import WCS

    # 1) Quick size check
    if image_a.shape != image_b.shape:
        return False

    # Extract WCS objects if headers provided
    if isinstance(header_a, WCS):
        wcs_a = header_a
        hdr_a = None
    else:
        hdr_a = header_a
        try:
            wcs_a = WCS(header_a)
        except Exception:
            return False

    if isinstance(header_b, WCS):
        wcs_b = header_b
        hdr_b = None
    else:
        hdr_b = header_b
        try:
            wcs_b = WCS(header_b)
        except Exception:
            return False

    # 2) Check equinox if headers available
    if hdr_a is not None and hdr_b is not None:
        equinox_a = get_equinox(hdr_a)
        equinox_b = get_equinox(hdr_b)
        if equinox_a is not None and equinox_b is not None:
            if equinox_a != equinox_b:
                return False

    # 2a) Check plate scale
    try:
        cdelt_a = get_platescale(wcs_a)
        cdelt_b = get_platescale(wcs_b)

        if abs(cdelt_a[0] - cdelt_b[0]) > tolerance:
            return False
        if abs(cdelt_a[1] - cdelt_b[1]) > tolerance:
            return False
    except ValueError:
        # Can't determine plate scale - continue with coordinate check
        pass

    # 2b) Compare coordinates of all pixels (brute-force for thorough check)
    # For large images, sample instead
    ny, nx = image_a.shape

    # Create pixel coordinate grids
    if nx * ny > 1000000:
        # Sample for large images
        x_sample = np.linspace(0, nx - 1, min(100, nx)).astype(int)
        y_sample = np.linspace(0, ny - 1, min(100, ny)).astype(int)
        x_grid, y_grid = np.meshgrid(x_sample, y_sample)
    else:
        y_grid, x_grid = np.mgrid[:ny, :nx]

    x_flat = x_grid.flatten()
    y_flat = y_grid.flatten()

    # Convert pixel to world coordinates
    try:
        world_a = wcs_a.pixel_to_world_values(x_flat, y_flat)
        world_b = wcs_b.pixel_to_world_values(x_flat, y_flat)

        ra_a, dec_a = world_a
        ra_b, dec_b = world_b

        # Check RA differences
        ra_diff = np.abs(ra_a - ra_b)
        # Handle RA wrap-around
        ra_diff = np.minimum(ra_diff, 360.0 - ra_diff)
        if np.any(ra_diff > tolerance):
            return False

        # Check Dec differences
        dec_diff = np.abs(dec_a - dec_b)
        if np.any(dec_diff > tolerance):
            return False

    except Exception:
        # WCS transformation failed
        return False

    return True


def images_aligned(
    image_a: np.ndarray,
    wcs_a: 'WCS',
    image_b: np.ndarray,
    wcs_b: 'WCS',
    tolerance_pixels: float = 0.5,
) -> bool:
    """
    Check if two images are aligned (same pixel grid).

    This is a more permissive check than astrometry_equal - it checks
    if the images can be directly compared pixel-by-pixel.

    Args:
        image_a: First image array
        wcs_a: WCS for first image
        image_b: Second image array
        wcs_b: WCS for second image
        tolerance_pixels: Allowable pixel offset

    Returns:
        True if images are aligned
    """
    # Same shape required
    if image_a.shape != image_b.shape:
        return False

    ny, nx = image_a.shape

    # Check corner pixels
    corners_x = [0, nx - 1, 0, nx - 1]
    corners_y = [0, 0, ny - 1, ny - 1]

    try:
        # Get world coordinates of corners from image A
        world_corners = wcs_a.pixel_to_world_values(corners_x, corners_y)

        # Convert back to pixel coordinates in image B
        pixels_b = wcs_b.world_to_pixel_values(*world_corners)

        # Check if pixel coordinates match
        x_b, y_b = pixels_b
        for i in range(4):
            if abs(x_b[i] - corners_x[i]) > tolerance_pixels:
                return False
            if abs(y_b[i] - corners_y[i]) > tolerance_pixels:
                return False

    except Exception:
        return False

    return True


def compute_pixel_offset(
    wcs_from: 'WCS',
    wcs_to: 'WCS',
    x: float,
    y: float,
) -> Tuple[float, float]:
    """
    Compute pixel offset between two WCS systems.

    Args:
        wcs_from: Source WCS
        wcs_to: Target WCS
        x: X pixel coordinate in source
        y: Y pixel coordinate in source

    Returns:
        Tuple of (dx, dy) pixel offset
    """
    # Convert pixel to world in source WCS
    world = wcs_from.pixel_to_world_values(x, y)

    # Convert world to pixel in target WCS
    x_to, y_to = wcs_to.world_to_pixel_values(*world)

    return x_to - x, y_to - y

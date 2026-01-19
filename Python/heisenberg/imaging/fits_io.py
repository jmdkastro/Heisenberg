"""
FITS file I/O utilities for Heisenberg.

This module provides functions for reading and writing FITS files,
extracting header information, and handling WCS astrometry.
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional, Dict, Any, Union
from dataclasses import dataclass

from astropy.io import fits
from astropy.wcs import WCS


@dataclass
class FitsImage:
    """Container for a FITS image with header metadata."""

    data: np.ndarray
    header: fits.Header
    filepath: Optional[Path] = None

    @property
    def wcs(self) -> WCS:
        """Get WCS object from header."""
        return WCS(self.header, naxis=2)

    @property
    def shape(self) -> Tuple[int, ...]:
        """Get image shape."""
        return self.data.shape

    @property
    def bmaj(self) -> Optional[float]:
        """Get beam major axis in degrees."""
        return self.header.get('BMAJ')

    @property
    def bmin(self) -> Optional[float]:
        """Get beam minor axis in degrees."""
        return self.header.get('BMIN')

    @property
    def beam(self) -> Optional[float]:
        """Get geometric mean beam size in degrees."""
        if self.bmaj is not None and self.bmin is not None:
            return np.sqrt(self.bmaj * self.bmin)
        return None

    @property
    def beam_arcsec(self) -> Optional[float]:
        """Get geometric mean beam size in arcseconds."""
        if self.beam is not None:
            return self.beam * 3600.0
        return None


def read_fits(filepath: Union[str, Path], extension: int = 0) -> FitsImage:
    """
    Read a FITS file and return data with header.

    Args:
        filepath: Path to the FITS file
        extension: HDU extension to read (default 0)

    Returns:
        FitsImage containing data array and header
    """
    filepath = Path(filepath)
    with fits.open(filepath) as hdul:
        data = hdul[extension].data.astype(np.float64)
        header = hdul[extension].header.copy()

    return FitsImage(data=data, header=header, filepath=filepath)


def write_fits(
    filepath: Union[str, Path],
    data: np.ndarray,
    header: Optional[fits.Header] = None,
    overwrite: bool = True
) -> None:
    """
    Write a FITS file.

    Args:
        filepath: Path to write the FITS file
        data: Image data array
        header: FITS header (optional)
        overwrite: Whether to overwrite existing file
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(filepath, overwrite=overwrite)


def get_platescale(header: fits.Header, tolerance: float = 1e-6) -> float:
    """
    Get the plate scale (pixel size) from a FITS header in degrees.

    Tries CDELT1/CDELT2 first, then CD matrix.

    Args:
        header: FITS header
        tolerance: Tolerance for comparing x and y scales

    Returns:
        Plate scale in degrees per pixel (absolute value)

    Raises:
        ValueError: If plate scale cannot be determined or x/y scales differ
    """
    # Try CDELT keywords first
    cdelt1 = header.get('CDELT1')
    cdelt2 = header.get('CDELT2')

    if cdelt1 is not None and cdelt2 is not None:
        cdelt1, cdelt2 = abs(cdelt1), abs(cdelt2)
        if abs(cdelt1 - cdelt2) > tolerance:
            raise ValueError(
                f"Non-square pixels: CDELT1={cdelt1}, CDELT2={cdelt2}"
            )
        return cdelt1

    # Try CD matrix
    cd1_1 = header.get('CD1_1')
    cd2_2 = header.get('CD2_2')
    cd1_2 = header.get('CD1_2', 0.0)
    cd2_1 = header.get('CD2_1', 0.0)

    if cd1_1 is not None and cd2_2 is not None:
        # For non-rotated images, CD1_1 and CD2_2 are the plate scales
        scale1 = np.sqrt(cd1_1**2 + cd2_1**2)
        scale2 = np.sqrt(cd1_2**2 + cd2_2**2)
        if abs(scale1 - scale2) > tolerance:
            raise ValueError(
                f"Non-square pixels from CD matrix: {scale1}, {scale2}"
            )
        return scale1

    raise ValueError("Cannot determine plate scale from header")


def get_platescale_signed(
    header: fits.Header,
    axis: str = 'x'
) -> float:
    """
    Get the signed plate scale for a specific axis.

    Args:
        header: FITS header
        axis: 'x' or 'y'

    Returns:
        Signed plate scale in degrees per pixel
    """
    if axis.lower() == 'x':
        cdelt = header.get('CDELT1')
        if cdelt is not None:
            return cdelt
        cd = header.get('CD1_1')
        if cd is not None:
            return cd
    elif axis.lower() == 'y':
        cdelt = header.get('CDELT2')
        if cdelt is not None:
            return cdelt
        cd = header.get('CD2_2')
        if cd is not None:
            return cd

    raise ValueError(f"Cannot determine {axis} plate scale from header")


def astrometry_equal(
    img1: FitsImage,
    img2: FitsImage,
    tolerance: float = 1e-6
) -> bool:
    """
    Check if two images have the same astrometry.

    Compares NAXIS, CRPIX, CRVAL, and CDELT/CD matrix values.

    Args:
        img1: First FITS image
        img2: Second FITS image
        tolerance: Tolerance for floating point comparisons (degrees)

    Returns:
        True if astrometry matches within tolerance
    """
    h1, h2 = img1.header, img2.header

    # Check dimensions
    if h1.get('NAXIS1') != h2.get('NAXIS1'):
        return False
    if h1.get('NAXIS2') != h2.get('NAXIS2'):
        return False

    # Check reference pixel
    if abs(h1.get('CRPIX1', 0) - h2.get('CRPIX1', 0)) > 0.5:
        return False
    if abs(h1.get('CRPIX2', 0) - h2.get('CRPIX2', 0)) > 0.5:
        return False

    # Check reference value
    if abs(h1.get('CRVAL1', 0) - h2.get('CRVAL1', 0)) > tolerance:
        return False
    if abs(h1.get('CRVAL2', 0) - h2.get('CRVAL2', 0)) > tolerance:
        return False

    # Check plate scale
    try:
        scale1 = get_platescale(h1, tolerance)
        scale2 = get_platescale(h2, tolerance)
        if abs(scale1 - scale2) > tolerance:
            return False
    except ValueError:
        return False

    return True


def update_header_beam(
    header: fits.Header,
    beam_size: float
) -> fits.Header:
    """
    Update header with new beam size.

    Args:
        header: Original FITS header
        beam_size: New beam size in arcseconds

    Returns:
        Updated header copy
    """
    new_header = header.copy()
    beam_deg = beam_size / 3600.0
    new_header['BMAJ'] = beam_deg
    new_header['BMIN'] = beam_deg
    return new_header


def update_header_regrid(
    header: fits.Header,
    scale_factor: float
) -> fits.Header:
    """
    Update header for regridded (resampled) image.

    Args:
        header: Original FITS header
        scale_factor: Factor by which pixel size increased

    Returns:
        Updated header copy
    """
    new_header = header.copy()

    # Update dimensions
    naxis1 = header.get('NAXIS1', 1)
    naxis2 = header.get('NAXIS2', 1)
    new_header['NAXIS1'] = int(round(naxis1 / scale_factor))
    new_header['NAXIS2'] = int(round(naxis2 / scale_factor))

    # Update reference pixel
    crpix1 = header.get('CRPIX1', 1)
    crpix2 = header.get('CRPIX2', 1)
    new_header['CRPIX1'] = round(crpix1 / scale_factor)
    new_header['CRPIX2'] = round(crpix2 / scale_factor)

    # Update plate scale
    if 'CDELT1' in header:
        new_header['CDELT1'] = header['CDELT1'] * scale_factor
    if 'CDELT2' in header:
        new_header['CDELT2'] = header['CDELT2'] * scale_factor
    if 'CD1_1' in header:
        new_header['CD1_1'] = header['CD1_1'] * scale_factor
    if 'CD2_2' in header:
        new_header['CD2_2'] = header['CD2_2'] * scale_factor

    return new_header

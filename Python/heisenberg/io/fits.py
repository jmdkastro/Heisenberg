"""
FITS file I/O for Heisenberg.

Provides functions for reading and writing FITS files with WCS support,
using astropy.

References:
    IDL: readfits, writefits, read_header
"""

from pathlib import Path
from typing import Optional, Tuple, Dict, Any, Union
import numpy as np

try:
    from astropy.io import fits
    from astropy.wcs import WCS
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False


def check_astropy():
    """Raise ImportError if astropy is not available."""
    if not HAS_ASTROPY:
        raise ImportError(
            "astropy required for FITS I/O: pip install astropy"
        )


def read_fits(
    path: Union[str, Path],
    extension: int = 0,
) -> Tuple[np.ndarray, Optional["WCS"], dict]:
    """
    Read a FITS file.

    Args:
        path: Path to FITS file
        extension: HDU extension to read (default: 0)

    Returns:
        Tuple of (data, wcs, header_dict)
        - data: Image array
        - wcs: WCS object if available, None otherwise
        - header_dict: Dictionary of header keywords
    """
    check_astropy()

    path = Path(path)
    with fits.open(path) as hdul:
        hdu = hdul[extension]
        data = hdu.data
        header = hdu.header

        # Try to extract WCS
        try:
            wcs = WCS(header)
            if wcs.naxis == 0:
                wcs = None
        except Exception:
            wcs = None

        # Convert header to dict
        header_dict = dict(header)

    return data, wcs, header_dict


def read_fits_data(
    path: Union[str, Path],
    extension: int = 0,
) -> np.ndarray:
    """
    Read only the data array from a FITS file.

    Args:
        path: Path to FITS file
        extension: HDU extension to read

    Returns:
        Image array
    """
    data, _, _ = read_fits(path, extension)
    return data


def read_fits_header(
    path: Union[str, Path],
    extension: int = 0,
) -> dict:
    """
    Read only the header from a FITS file.

    Args:
        path: Path to FITS file
        extension: HDU extension to read

    Returns:
        Dictionary of header keywords
    """
    check_astropy()

    path = Path(path)
    with fits.open(path) as hdul:
        header = dict(hdul[extension].header)

    return header


def write_fits(
    data: np.ndarray,
    path: Union[str, Path],
    header: Optional[Dict[str, Any]] = None,
    wcs: Optional["WCS"] = None,
    overwrite: bool = False,
) -> None:
    """
    Write data to a FITS file.

    Args:
        data: Image array to write
        path: Output file path
        header: Optional header keywords
        wcs: Optional WCS object
        overwrite: Overwrite existing file
    """
    check_astropy()

    path = Path(path)

    # Create header
    hdu_header = fits.Header()

    # Add WCS if provided
    if wcs is not None:
        hdu_header.update(wcs.to_header())

    # Add custom header keywords
    if header is not None:
        for key, value in header.items():
            # Skip reserved keywords
            if key.upper() in ('SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2'):
                continue
            try:
                hdu_header[key] = value
            except ValueError:
                # Skip invalid keywords
                pass

    hdu = fits.PrimaryHDU(data, header=hdu_header)
    hdu.writeto(path, overwrite=overwrite)


def get_pixel_scale(
    wcs: "WCS",
) -> float:
    """
    Get pixel scale from WCS in degrees per pixel.

    Args:
        wcs: WCS object

    Returns:
        Pixel scale in degrees per pixel
    """
    check_astropy()

    if wcs is None:
        raise ValueError("No WCS provided")

    # Get pixel scale from CD matrix or CDELT
    try:
        # Try CD matrix first
        if hasattr(wcs.wcs, 'cd') and wcs.wcs.cd is not None:
            cd = wcs.wcs.cd
            scale = np.sqrt(np.abs(cd[0, 0] * cd[1, 1] - cd[0, 1] * cd[1, 0]))
        else:
            # Fall back to CDELT
            scale = np.abs(wcs.wcs.cdelt[0])
    except Exception:
        raise ValueError("Cannot determine pixel scale from WCS")

    return scale


def get_pixel_scale_pc(
    wcs: "WCS",
    distance_pc: float,
) -> float:
    """
    Get pixel scale in parsecs per pixel.

    Args:
        wcs: WCS object
        distance_pc: Distance to source in parsecs

    Returns:
        Pixel scale in parsecs per pixel
    """
    scale_deg = get_pixel_scale(wcs)
    scale_rad = np.deg2rad(scale_deg)
    scale_pc = scale_rad * distance_pc

    return scale_pc


def load_tracer_maps(
    star_path: Union[str, Path],
    gas_path: Union[str, Path],
    distance_pc: float,
    mask_path: Optional[Union[str, Path]] = None,
) -> Tuple[np.ndarray, np.ndarray, float, Optional[np.ndarray]]:
    """
    Load stellar and gas tracer maps.

    Args:
        star_path: Path to stellar tracer FITS file
        gas_path: Path to gas tracer FITS file
        distance_pc: Distance to galaxy in parsecs
        mask_path: Optional path to mask FITS file

    Returns:
        Tuple of (star_image, gas_image, pixel_scale_pc, mask)
    """
    star_data, star_wcs, _ = read_fits(star_path)
    gas_data, gas_wcs, _ = read_fits(gas_path)

    # Check shape compatibility
    if star_data.shape != gas_data.shape:
        raise ValueError(
            f"Image shape mismatch: star={star_data.shape}, gas={gas_data.shape}"
        )

    # Get pixel scale
    if star_wcs is not None:
        pixel_scale = get_pixel_scale_pc(star_wcs, distance_pc)
    elif gas_wcs is not None:
        pixel_scale = get_pixel_scale_pc(gas_wcs, distance_pc)
    else:
        raise ValueError("No WCS found in either file - cannot determine pixel scale")

    # Load mask if provided
    mask = None
    if mask_path is not None:
        mask, _, _ = read_fits(mask_path)
        if mask.shape != star_data.shape:
            raise ValueError(
                f"Mask shape mismatch: mask={mask.shape}, data={star_data.shape}"
            )

    return star_data, gas_data, pixel_scale, mask


def load_sensitivity_maps(
    star_sens_path: Optional[Union[str, Path]] = None,
    gas_sens_path: Optional[Union[str, Path]] = None,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Load sensitivity maps for star and gas tracers.

    Args:
        star_sens_path: Path to stellar sensitivity FITS file
        gas_sens_path: Path to gas sensitivity FITS file

    Returns:
        Tuple of (star_sensitivity, gas_sensitivity)
    """
    star_sens = None
    gas_sens = None

    if star_sens_path is not None:
        star_sens, _, _ = read_fits(star_sens_path)

    if gas_sens_path is not None:
        gas_sens, _, _ = read_fits(gas_sens_path)

    return star_sens, gas_sens

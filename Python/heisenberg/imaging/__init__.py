"""
Imaging module for Heisenberg.

Contains FITS I/O, image smoothing/convolution, masking utilities,
and aperture photometry functions.
"""

from heisenberg.imaging.fits_io import (
    FitsImage,
    read_fits,
    write_fits,
    get_platescale,
    get_platescale_signed,
    astrometry_equal,
    update_header_beam,
    update_header_regrid,
)
from heisenberg.imaging.convolution import (
    psf_gaussian,
    psf_tophat,
    convolve_image,
    smooth_to_resolution,
    regrid_image,
)
from heisenberg.imaging.masking import (
    Mask,
    create_circular_mask,
    create_elliptical_mask,
    create_radial_mask,
    create_box_mask,
    create_polygon_mask,
    mask_from_nan,
    synchronize_masks,
    apply_mask_to_image,
    write_mask_fits,
)
from heisenberg.imaging.aperture import (
    ApertureFluxResult,
    generate_aperture_sizes,
    create_smoothed_cube,
    measure_aperture_flux,
    compute_flux_ratios,
    select_non_overlapping_peaks,
    monte_carlo_flux_ratios,
)
from heisenberg.imaging.mask_tool import (
    MaskToolResult,
    mask_tool,
    create_mask_from_ds9,
    region_to_mask,
    regions_to_mask,
)
from heisenberg.imaging.astrometry import (
    get_equinox,
    get_rotation,
    astrometry_equal as astrometry_equal_thorough,
    images_aligned,
    compute_pixel_offset,
)

__all__ = [
    # FITS I/O
    "FitsImage",
    "read_fits",
    "write_fits",
    "get_platescale",
    "get_platescale_signed",
    "astrometry_equal",
    "update_header_beam",
    "update_header_regrid",
    # Convolution
    "psf_gaussian",
    "psf_tophat",
    "convolve_image",
    "smooth_to_resolution",
    "regrid_image",
    # Masking
    "Mask",
    "create_circular_mask",
    "create_elliptical_mask",
    "create_radial_mask",
    "create_box_mask",
    "create_polygon_mask",
    "mask_from_nan",
    "synchronize_masks",
    "apply_mask_to_image",
    "write_mask_fits",
    # Aperture photometry
    "ApertureFluxResult",
    "generate_aperture_sizes",
    "create_smoothed_cube",
    "measure_aperture_flux",
    "compute_flux_ratios",
    "select_non_overlapping_peaks",
    "monte_carlo_flux_ratios",
    # Mask tool
    "MaskToolResult",
    "mask_tool",
    "create_mask_from_ds9",
    "region_to_mask",
    "regions_to_mask",
    # Astrometry
    "get_equinox",
    "get_rotation",
    "astrometry_equal_thorough",
    "images_aligned",
    "compute_pixel_offset",
]

"""
High-level masking tool for Heisenberg.

Provides unified workflow for applying DS9 region masks to astronomical
images. This is a Python translation of IDL mask_tool.pro.

Usage:
    >>> from heisenberg.imaging.mask_tool import mask_tool
    >>> masked_image, mask = mask_tool(
    ...     image_path='galaxy.fits',
    ...     ds9_positive_path='allowed_regions.reg',
    ...     ds9_negative_path='blocked_regions.reg',
    ... )

References:
    IDL: mask_tool.pro, mask_ds9_file_mask.pro
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np

# Import directly from modules to avoid circular imports
from heisenberg.imaging.masking import (
    Mask,
    create_circular_mask,
    create_elliptical_mask,
    create_polygon_mask,
    write_mask_fits,
)


@dataclass
class MaskToolResult:
    """
    Result from mask_tool operation.

    Attributes:
        masked_image: Image with masked pixels set to NaN
        mask: Mask array (1.0 = unmasked, 0.0 = masked)
        header: FITS header from the input image
    """
    masked_image: np.ndarray
    mask: Mask
    header: Optional[dict] = None


def region_to_mask(
    region: 'DS9Region',
    shape: Tuple[int, int],
    negative: bool = False,
) -> Mask:
    """
    Convert a DS9Region to a Mask.

    Args:
        region: DS9Region object
        shape: Image shape (ny, nx)
        negative: If True, mask inside region; if False, keep inside region

    Returns:
        Mask object
    """
    # Import here to avoid circular import
    from heisenberg.plotting.ds9_regions import RegionType, box_vertices

    ny, nx = shape
    params = region.params

    # DS9 uses 1-indexed coordinates, convert to 0-indexed
    # Also swap x/y for numpy array indexing (row, col)

    if region.region_type == RegionType.CIRCLE:
        x = params['x'] - 1  # Convert to 0-indexed
        y = params['y'] - 1
        radius = params['radius']

        # create_circular_mask expects (y, x) center for numpy indexing
        mask = create_circular_mask(
            shape=shape,
            center=(y, x),
            radius=radius,
            inside=negative,  # If negative, mask inside
        )

    elif region.region_type == RegionType.ELLIPSE:
        x = params['x'] - 1
        y = params['y'] - 1
        semi_major = params['semi_major']
        semi_minor = params['semi_minor']
        angle_deg = params['angle']

        # Convert angle to radians (DS9 uses degrees)
        angle_rad = np.deg2rad(angle_deg)

        mask = create_elliptical_mask(
            shape=shape,
            center=(y, x),
            semi_major=semi_major,
            semi_minor=semi_minor,
            position_angle=angle_rad,
            inside=negative,
        )

    elif region.region_type == RegionType.BOX:
        x = params['x'] - 1
        y = params['y'] - 1
        width = params['width']
        height = params['height']
        angle_deg = params['angle']

        # Get box vertices
        angle_rad = np.deg2rad(angle_deg)
        x_verts, y_verts = box_vertices(x, y, width, height, angle_rad)

        # Convert to (row, col) corner coordinates for polygon mask
        corners = [(yv, xv) for xv, yv in zip(x_verts, y_verts)]

        mask = create_polygon_mask(
            shape=shape,
            vertices=corners,
            inside=negative,
        )

    elif region.region_type == RegionType.POLYGON:
        x_coords = [xc - 1 for xc in params['x']]  # Convert to 0-indexed
        y_coords = [yc - 1 for yc in params['y']]

        # Convert to (row, col) vertices
        vertices = [(y, x) for x, y in zip(x_coords, y_coords)]

        mask = create_polygon_mask(
            shape=shape,
            vertices=vertices,
            inside=negative,
        )

    else:
        # Unsupported region type - return all-unmasked
        mask = Mask(
            data=np.ones(shape),
            description=f"Unsupported region type: {region.region_type}",
        )

    return mask


def regions_to_mask(
    regions: list,
    shape: Tuple[int, int],
    negative: bool = False,
) -> Mask:
    """
    Convert multiple DS9 regions to a combined mask.

    For positive masks (negative=False): regions define allowed areas.
    For negative masks (negative=True): regions define blocked areas.

    Args:
        regions: List of DS9Region objects
        shape: Image shape (ny, nx)
        negative: If True, regions are blocked; if False, regions are allowed

    Returns:
        Combined Mask object
    """
    if not regions:
        # No regions - return appropriate default
        if negative:
            return Mask(data=np.ones(shape), description="No negative regions")
        else:
            return Mask(data=np.zeros(shape), description="No positive regions")

    # Start with appropriate base
    if negative:
        # Negative mask: start unmasked, then block each region
        combined = np.ones(shape)
        for region in regions:
            if region.is_exclude:
                # Excluded region in negative mask = allow this region
                region_mask = region_to_mask(region, shape, negative=True)
                combined = np.maximum(combined, 1.0 - region_mask.data)
            else:
                # Normal region in negative mask = block this region
                region_mask = region_to_mask(region, shape, negative=True)
                combined = np.minimum(combined, region_mask.data)
    else:
        # Positive mask: start masked, then allow each region
        combined = np.zeros(shape)
        for region in regions:
            if region.is_exclude:
                # Excluded region in positive mask = block this region
                pass  # Handled after combining positives
            else:
                # Normal region in positive mask = allow this region
                # region_to_mask with negative=False returns mask where inside=1 (unmasked)
                region_mask = region_to_mask(region, shape, negative=False)
                combined = np.maximum(combined, region_mask.data)

    desc = f"{'Negative' if negative else 'Positive'} mask from {len(regions)} regions"
    return Mask(data=combined, description=desc)


def mask_tool(
    image_input: Union[str, Path, np.ndarray],
    image_header: Optional[dict] = None,
    ds9_positive_path: Optional[Union[str, Path]] = None,
    ds9_negative_path: Optional[Union[str, Path]] = None,
    masked_image_path: Optional[Union[str, Path]] = None,
    mask_path: Optional[Union[str, Path]] = None,
    convert: bool = False,
    conv_filepath: Optional[Union[str, Path]] = None,
    run_without_masks: bool = False,
) -> MaskToolResult:
    """
    Apply DS9 region masks to an astronomical image.

    This is a Python translation of IDL mask_tool.pro.

    Args:
        image_input: Path to FITS file or image array
        image_header: FITS header (required if image_input is array)
        ds9_positive_path: Path to DS9 region file for allowed regions
        ds9_negative_path: Path to DS9 region file for blocked regions
        masked_image_path: Path to output masked image FITS file
        mask_path: Path to output mask FITS file
        convert: If True, convert DS9 files to image coordinates first
        conv_filepath: Path to FITS file for coordinate conversion
        run_without_masks: If True, allow running without any mask files

    Returns:
        MaskToolResult with masked image, mask array, and header

    Raises:
        ValueError: If no mask files provided and run_without_masks=False
        FileNotFoundError: If DS9 conversion requested but DS9 not installed

    Example:
        >>> result = mask_tool(
        ...     'galaxy.fits',
        ...     ds9_positive_path='allowed.reg',
        ...     ds9_negative_path='blocked.reg',
        ... )
        >>> masked = result.masked_image
        >>> mask = result.mask
    """
    from heisenberg.io.fits import read_fits, write_fits

    # Load image
    if isinstance(image_input, (str, Path)):
        image_path = Path(image_input)
        image, _, header = read_fits(image_path)
        if conv_filepath is None:
            conv_filepath = image_path
    else:
        image = image_input
        header = image_header
        if convert and conv_filepath is None:
            raise ValueError(
                "conv_filepath required for coordinate conversion "
                "when image_input is an array"
            )

    shape = image.shape

    # Check that at least one mask is provided
    if ds9_positive_path is None and ds9_negative_path is None:
        if run_without_masks:
            # No masks - return image unchanged with all-ones mask
            mask = Mask(data=np.ones(shape), description="No masks applied")
            return MaskToolResult(
                masked_image=image.copy(),
                mask=mask,
                header=header,
            )
        else:
            raise ValueError(
                "Neither positive nor negative DS9 mask provided. "
                "Supply at least one or set run_without_masks=True."
            )

    # Convert DS9 files to image coordinates if requested
    # Import here to avoid circular import
    from heisenberg.plotting.ds9_regions import (
        parse_ds9_region_file,
        ds9_convert_to_image,
    )

    if convert:
        if ds9_positive_path is not None:
            ds9_positive_path = Path(ds9_positive_path)
            positive_converted = ds9_positive_path.with_suffix('.converted.reg')
            success = ds9_convert_to_image(
                conv_filepath, ds9_positive_path, positive_converted
            )
            if success:
                ds9_positive_path = positive_converted

        if ds9_negative_path is not None:
            ds9_negative_path = Path(ds9_negative_path)
            negative_converted = ds9_negative_path.with_suffix('.converted.reg')
            success = ds9_convert_to_image(
                conv_filepath, ds9_negative_path, negative_converted
            )
            if success:
                ds9_negative_path = negative_converted

    # Parse and create masks
    positive_mask = None
    negative_mask = None

    if ds9_positive_path is not None:
        regions, _ = parse_ds9_region_file(ds9_positive_path)
        positive_mask = regions_to_mask(regions, shape, negative=False)

    if ds9_negative_path is not None:
        regions, _ = parse_ds9_region_file(ds9_negative_path)
        negative_mask = regions_to_mask(regions, shape, negative=True)

    # Combine masks
    if positive_mask is not None and negative_mask is not None:
        # Both masks: only regions allowed by BOTH are unmasked
        combined_data = positive_mask.data * negative_mask.data
        mask = Mask(
            data=combined_data,
            description="Combined positive and negative masks",
        )
    elif positive_mask is not None:
        mask = positive_mask
    else:
        mask = negative_mask

    # Apply mask to image
    masked_image = mask.apply(image)

    # Write outputs if paths provided
    if mask_path is not None:
        write_mask_fits(mask_path, mask, header=header)

    if masked_image_path is not None:
        write_fits(masked_image, masked_image_path, header=header, overwrite=True)

    return MaskToolResult(
        masked_image=masked_image,
        mask=mask,
        header=header,
    )


def create_mask_from_ds9(
    ds9_path: Union[str, Path],
    shape: Tuple[int, int],
    negative: bool = False,
) -> Mask:
    """
    Create a mask from a DS9 region file.

    Convenience function for creating masks without the full mask_tool workflow.

    Args:
        ds9_path: Path to DS9 region file
        shape: Image shape (ny, nx)
        negative: If True, regions are blocked; if False, regions are allowed

    Returns:
        Mask object
    """
    # Import here to avoid circular import
    from heisenberg.plotting.ds9_regions import parse_ds9_region_file

    regions, _ = parse_ds9_region_file(ds9_path)
    return regions_to_mask(regions, shape, negative=negative)

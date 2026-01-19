"""
Image masking utilities for Heisenberg.

This module provides functions for creating and applying masks to
astronomical images, including support for DS9 region files and
radial cuts.
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional, Union, List
from dataclasses import dataclass

from astropy.io import fits


@dataclass
class Mask:
    """Container for an image mask."""

    data: np.ndarray  # 1.0 = unmasked, 0.0 = masked
    description: str = ""

    @property
    def shape(self) -> Tuple[int, ...]:
        """Get mask shape."""
        return self.data.shape

    @property
    def n_unmasked(self) -> int:
        """Number of unmasked pixels."""
        return int(np.sum(self.data > 0.5))

    @property
    def n_masked(self) -> int:
        """Number of masked pixels."""
        return int(np.sum(self.data < 0.5))

    @property
    def fraction_unmasked(self) -> float:
        """Fraction of unmasked pixels."""
        return self.n_unmasked / self.data.size

    def apply(self, image: np.ndarray) -> np.ndarray:
        """
        Apply mask to image, setting masked pixels to NaN.

        Args:
            image: Input image

        Returns:
            Image with masked pixels set to NaN
        """
        result = image.copy()
        result[self.data < 0.5] = np.nan
        return result

    def combine(self, other: 'Mask') -> 'Mask':
        """
        Combine with another mask (logical AND).

        Args:
            other: Another mask of the same shape

        Returns:
            Combined mask where both masks must be unmasked
        """
        if self.shape != other.shape:
            raise ValueError("Masks must have the same shape")
        return Mask(
            data=self.data * other.data,
            description=f"({self.description}) AND ({other.description})"
        )


def create_circular_mask(
    shape: Tuple[int, int],
    center: Tuple[float, float],
    radius: float,
    inside: bool = True
) -> Mask:
    """
    Create a circular mask.

    Args:
        shape: Image shape (ny, nx)
        center: Center coordinates (y, x) in pixels
        radius: Radius in pixels
        inside: If True, mask inside circle; if False, mask outside

    Returns:
        Mask with specified region masked
    """
    ny, nx = shape
    y, x = np.ogrid[:ny, :nx]
    cy, cx = center

    dist = np.sqrt((x - cx)**2 + (y - cy)**2)

    if inside:
        mask_data = np.where(dist <= radius, 0.0, 1.0)
        desc = f"Circle r={radius:.1f} masked inside"
    else:
        mask_data = np.where(dist <= radius, 1.0, 0.0)
        desc = f"Circle r={radius:.1f} masked outside"

    return Mask(data=mask_data, description=desc)


def create_elliptical_mask(
    shape: Tuple[int, int],
    center: Tuple[float, float],
    semi_major: float,
    semi_minor: float,
    position_angle: float = 0.0,
    inside: bool = True
) -> Mask:
    """
    Create an elliptical mask.

    Args:
        shape: Image shape (ny, nx)
        center: Center coordinates (y, x) in pixels
        semi_major: Semi-major axis in pixels
        semi_minor: Semi-minor axis in pixels
        position_angle: Position angle in radians (from y-axis)
        inside: If True, mask inside ellipse; if False, mask outside

    Returns:
        Mask with specified region masked
    """
    ny, nx = shape
    y, x = np.ogrid[:ny, :nx]
    cy, cx = center

    # Rotate coordinates
    cos_pa = np.cos(position_angle)
    sin_pa = np.sin(position_angle)
    x_rot = (x - cx) * cos_pa + (y - cy) * sin_pa
    y_rot = -(x - cx) * sin_pa + (y - cy) * cos_pa

    # Ellipse equation
    ellipse = (x_rot / semi_major)**2 + (y_rot / semi_minor)**2

    if inside:
        mask_data = np.where(ellipse <= 1.0, 0.0, 1.0)
        desc = f"Ellipse a={semi_major:.1f} b={semi_minor:.1f} masked inside"
    else:
        mask_data = np.where(ellipse <= 1.0, 1.0, 0.0)
        desc = f"Ellipse a={semi_major:.1f} b={semi_minor:.1f} masked outside"

    return Mask(data=mask_data, description=desc)


def create_radial_mask(
    shape: Tuple[int, int],
    center: Tuple[float, float],
    inner_radius: float,
    outer_radius: float,
    inclination: float = 0.0,
    position_angle: float = 0.0,
    pixel_to_pc: float = 1.0
) -> Mask:
    """
    Create a radial annulus mask with optional inclination correction.

    Args:
        shape: Image shape (ny, nx)
        center: Center coordinates (y, x) in pixels
        inner_radius: Inner radius in pc
        outer_radius: Outer radius in pc
        inclination: Inclination angle in radians
        position_angle: Position angle in radians
        pixel_to_pc: Conversion from pixels to pc

    Returns:
        Mask with pixels outside the radial range masked
    """
    ny, nx = shape
    y_idx, x_idx = np.ogrid[:ny, :nx]
    cy, cx = center

    # Distances in pixel coordinates
    dx_pix = x_idx - cx
    dy_pix = y_idx - cy

    # Rotate by position angle
    cos_pa = np.cos(-position_angle)
    sin_pa = np.sin(-position_angle)
    dx_rot = cos_pa * dx_pix - sin_pa * dy_pix
    dy_rot = sin_pa * dx_pix + cos_pa * dy_pix

    # Convert to physical coordinates with inclination correction
    dx_pc = dx_rot * pixel_to_pc
    dy_pc = dy_rot * pixel_to_pc / np.cos(inclination)

    # Calculate radius
    radius_pc = np.sqrt(dx_pc**2 + dy_pc**2)

    # Create mask
    mask_data = np.where(
        (radius_pc >= inner_radius) & (radius_pc <= outer_radius),
        1.0, 0.0
    )

    return Mask(
        data=mask_data,
        description=f"Radial mask {inner_radius:.0f}-{outer_radius:.0f} pc"
    )


def create_box_mask(
    shape: Tuple[int, int],
    corners: List[Tuple[float, float]],
    inside: bool = True
) -> Mask:
    """
    Create a rectangular/box mask.

    Args:
        shape: Image shape (ny, nx)
        corners: List of (y, x) corner coordinates
        inside: If True, mask inside box; if False, mask outside

    Returns:
        Mask with specified region masked
    """
    from matplotlib.path import Path as MplPath

    ny, nx = shape
    y, x = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    points = np.column_stack([y.ravel(), x.ravel()])

    # Create polygon path
    path = MplPath(corners)
    inside_mask = path.contains_points(points).reshape(shape)

    if inside:
        mask_data = np.where(inside_mask, 0.0, 1.0)
        desc = "Box masked inside"
    else:
        mask_data = np.where(inside_mask, 1.0, 0.0)
        desc = "Box masked outside"

    return Mask(data=mask_data, description=desc)


def create_polygon_mask(
    shape: Tuple[int, int],
    vertices: List[Tuple[float, float]],
    inside: bool = True
) -> Mask:
    """
    Create a polygon mask.

    Args:
        shape: Image shape (ny, nx)
        vertices: List of (y, x) vertex coordinates
        inside: If True, mask inside polygon; if False, mask outside

    Returns:
        Mask with specified region masked
    """
    from matplotlib.path import Path as MplPath

    ny, nx = shape
    y, x = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')
    points = np.column_stack([y.ravel(), x.ravel()])

    # Create polygon path
    path = MplPath(vertices)
    inside_mask = path.contains_points(points).reshape(shape)

    if inside:
        mask_data = np.where(inside_mask, 0.0, 1.0)
        desc = f"Polygon ({len(vertices)} vertices) masked inside"
    else:
        mask_data = np.where(inside_mask, 1.0, 0.0)
        desc = f"Polygon ({len(vertices)} vertices) masked outside"

    return Mask(data=mask_data, description=desc)


def mask_from_nan(image: np.ndarray) -> Mask:
    """
    Create a mask from NaN values in an image.

    Args:
        image: Input image

    Returns:
        Mask where NaN pixels are masked
    """
    mask_data = np.where(np.isnan(image), 0.0, 1.0)
    return Mask(data=mask_data, description="NaN mask")


def synchronize_masks(images: List[np.ndarray]) -> Mask:
    """
    Create a combined mask from NaN values in multiple images.

    A pixel is masked if it is NaN in ANY of the input images.

    Args:
        images: List of input images (must have same shape)

    Returns:
        Combined mask
    """
    if not images:
        raise ValueError("At least one image required")

    shape = images[0].shape
    for i, img in enumerate(images[1:], 2):
        if img.shape != shape:
            raise ValueError(f"Image {i} has different shape")

    # Start with all unmasked
    mask_data = np.ones(shape)

    # Mask where any image has NaN
    for img in images:
        mask_data[np.isnan(img)] = 0.0

    return Mask(
        data=mask_data,
        description=f"Synchronized mask from {len(images)} images"
    )


def apply_mask_to_image(
    image: np.ndarray,
    mask: Union[Mask, np.ndarray]
) -> np.ndarray:
    """
    Apply a mask to an image.

    Args:
        image: Input image
        mask: Mask object or array (1.0 = unmasked, 0.0 = masked)

    Returns:
        Image with masked pixels set to NaN
    """
    if isinstance(mask, Mask):
        return mask.apply(image)

    result = image.copy()
    result[mask < 0.5] = np.nan
    return result


def write_mask_fits(
    filepath: Union[str, Path],
    mask: Mask,
    header: Optional[fits.Header] = None,
    overwrite: bool = True
) -> None:
    """
    Write a mask to a FITS file.

    Args:
        filepath: Output file path
        mask: Mask to write
        header: FITS header (optional)
        overwrite: Whether to overwrite existing file
    """
    from heisenberg.imaging.fits_io import write_fits
    write_fits(filepath, mask.data, header=header, overwrite=overwrite)

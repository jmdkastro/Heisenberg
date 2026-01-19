"""
FITS file display functions.

Translates IDL plotfits.pro and plotfits_files.pro to Python.
"""

import math
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from heisenberg.plotting.image_display import display_image, get_percentile_scaling


def plot_fits(
    data: np.ndarray,
    log_scale: bool = True,
    log_offset: float = 0.0,
    right_offset: float = 0.0,
    vmin_percentile: float = 1.0,
    vmax_percentile: float = 99.0,
    cmap: str = 'viridis',
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (8, 6),
    colorbar: bool = True,
) -> plt.Axes:
    """
    Display a FITS image with automatic scaling.

    This function provides features similar to IDL's plotfits.pro:
    - Histogram-based contrast stretching via percentiles
    - Log/linear color scaling
    - Configurable offsets for min/max adjustment

    Args:
        data: 2D numpy array (image data)
        log_scale: If True, display log10(data)
        log_offset: Offset to add to log10(min), effectively raising floor
        right_offset: Offset to subtract from log10(max), lowering ceiling
        vmin_percentile: Lower percentile for auto scaling
        vmax_percentile: Upper percentile for auto scaling
        cmap: Matplotlib colormap name
        title: Plot title
        ax: Existing matplotlib axes (creates new if None)
        figsize: Figure size if creating new figure
        colorbar: If True, add colorbar

    Returns:
        The matplotlib Axes object
    """
    # Make a copy
    image = data.copy().astype(float)

    if log_scale:
        # Handle non-positive values for log
        positive_mask = image > 0
        if positive_mask.any():
            min_positive = np.nanmin(image[positive_mask])
        else:
            min_positive = 1e-10
        image[~positive_mask] = min_positive
        image[~np.isfinite(image)] = min_positive

        # Apply log transform
        log_image = np.log10(image)

        # Get scaling range with offsets
        vmin, vmax = get_percentile_scaling(
            log_image, vmin_percentile, vmax_percentile
        )
        vmin = vmin + log_offset
        vmax = vmax - right_offset

        # Display the log-transformed image using astropy stretch
        ax = display_image(
            log_image,
            stretch='linear',  # Already took log manually
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            title=title,
            ax=ax,
            figsize=figsize,
            colorbar=colorbar,
        )
    else:
        # Linear scaling
        vmin, vmax = get_percentile_scaling(
            image, vmin_percentile, vmax_percentile
        )
        # Apply offsets in log space then convert back
        if log_offset != 0 or right_offset != 0:
            log_vmin = np.log10(max(vmin, 1e-10)) + log_offset
            log_vmax = np.log10(max(vmax, 1e-10)) - right_offset
            vmin = 10**log_vmin
            vmax = 10**log_vmax

        ax = display_image(
            image,
            stretch='linear',
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            title=title,
            ax=ax,
            figsize=figsize,
            colorbar=colorbar,
        )

    return ax


def plot_fits_from_file(
    filepath: Union[str, Path],
    extension: int = 0,
    **kwargs
) -> plt.Axes:
    """
    Load and display a FITS file.

    Args:
        filepath: Path to FITS file
        extension: FITS extension number to read
        **kwargs: Additional arguments passed to plot_fits

    Returns:
        The matplotlib Axes object
    """
    from astropy.io import fits

    filepath = Path(filepath)

    # Handle files without .fits extension
    if not filepath.suffix:
        filepath = filepath.with_suffix('.fits')

    with fits.open(filepath) as hdul:
        data = hdul[extension].data

    # Use filename as default title
    if 'title' not in kwargs:
        kwargs['title'] = filepath.stem

    return plot_fits(data, **kwargs)


def plot_fits_grid(
    images: List[np.ndarray],
    ncols: int = 2,
    log_scale: bool = True,
    log_offsets: Optional[List[float]] = None,
    right_offsets: Optional[List[float]] = None,
    cmaps: Optional[List[str]] = None,
    titles: Optional[List[str]] = None,
    figsize: Optional[Tuple[float, float]] = None,
    shared_colorbar: bool = False,
) -> Figure:
    """
    Create multi-panel display of multiple images.

    This function provides features similar to IDL's plotfits_files.pro:
    - Grid layout with configurable columns
    - Per-image log offset and colormap
    - Optional shared or individual colorbars

    Args:
        images: List of 2D numpy arrays
        ncols: Number of columns in grid
        log_scale: If True, display log10(data)
        log_offsets: Per-image log offset (list or single value)
        right_offsets: Per-image right offset for max
        cmaps: Per-image colormap names
        titles: Per-image titles
        figsize: Figure size (auto-calculated if None)
        shared_colorbar: If True, use single colorbar for all panels

    Returns:
        The matplotlib Figure object
    """
    n_images = len(images)
    nrows = math.ceil(n_images / ncols)

    # Set defaults
    if log_offsets is None:
        log_offsets = [0.0] * n_images
    elif isinstance(log_offsets, (int, float)):
        log_offsets = [log_offsets] * n_images

    if right_offsets is None:
        right_offsets = [0.0] * n_images
    elif isinstance(right_offsets, (int, float)):
        right_offsets = [right_offsets] * n_images

    if cmaps is None:
        cmaps = ['viridis'] * n_images
    elif isinstance(cmaps, str):
        cmaps = [cmaps] * n_images

    if titles is None:
        titles = [f'Image {i+1}' for i in range(n_images)]

    # Calculate figure size
    if figsize is None:
        figsize = (4 * ncols, 3.5 * nrows)

    # Create figure and axes
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)

    # Flatten axes array for easy iteration
    if n_images == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Plot each image
    for i, (image, loff, roff, cmap, title) in enumerate(
        zip(images, log_offsets, right_offsets, cmaps, titles)
    ):
        ax = axes[i]
        plot_fits(
            image,
            log_scale=log_scale,
            log_offset=loff,
            right_offset=roff,
            cmap=cmap,
            title=title,
            ax=ax,
            colorbar=not shared_colorbar,
        )

    # Hide unused axes
    for i in range(n_images, len(axes)):
        axes[i].set_visible(False)

    # Add shared colorbar if requested
    if shared_colorbar and n_images > 0:
        # Use the last image's normalization for colorbar
        fig.colorbar(
            axes[0].images[0],
            ax=axes.tolist(),
            orientation='vertical',
            fraction=0.02,
            pad=0.04,
        )

    plt.tight_layout()
    return fig


def plot_fits_grid_from_files(
    filepaths: List[Union[str, Path]],
    extension: int = 0,
    **kwargs
) -> Figure:
    """
    Load and display multiple FITS files in a grid.

    Args:
        filepaths: List of paths to FITS files
        extension: FITS extension number to read
        **kwargs: Additional arguments passed to plot_fits_grid

    Returns:
        The matplotlib Figure object
    """
    from astropy.io import fits

    images = []
    for filepath in filepaths:
        filepath = Path(filepath)
        if not filepath.suffix:
            filepath = filepath.with_suffix('.fits')

        with fits.open(filepath) as hdul:
            images.append(hdul[extension].data.copy())

    # Use filenames as default titles
    if 'titles' not in kwargs:
        kwargs['titles'] = [Path(p).stem for p in filepaths]

    return plot_fits_grid(images, **kwargs)


def plot_star_gas_maps(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    star_title: str = 'Star Tracer',
    gas_title: str = 'Gas Tracer',
    log_scale: bool = True,
    figsize: Tuple[float, float] = (12, 5),
) -> Figure:
    """
    Display star and gas tracer maps side by side.

    This is a convenience function for the common use case of
    comparing star and gas images.

    Args:
        star_image: Stellar tracer map
        gas_image: Gas tracer map
        star_title: Title for star panel
        gas_title: Title for gas panel
        log_scale: If True, display in log scale
        figsize: Figure size

    Returns:
        The matplotlib Figure object
    """
    return plot_fits_grid(
        images=[star_image, gas_image],
        ncols=2,
        log_scale=log_scale,
        titles=[star_title, gas_title],
        figsize=figsize,
    )

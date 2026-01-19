"""
Core astronomical image display functions.

Thin wrappers around astropy.visualization for Heisenberg-specific needs.
Translates IDL disp.pro functionality.
"""

from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

try:
    from astropy.visualization import (
        simple_norm,
        imshow_norm,
        ZScaleInterval,
        MinMaxInterval,
        PercentileInterval,
        AsinhStretch,
        LogStretch,
        LinearStretch,
        SqrtStretch,
    )
    from astropy.wcs import WCS
    HAS_ASTROPY_VIS = True
except ImportError:
    HAS_ASTROPY_VIS = False
    WCS = None


def display_image(
    image: np.ndarray,
    wcs: Optional['WCS'] = None,
    stretch: str = 'linear',
    interval: str = 'minmax',
    percent: float = 99.0,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap: str = 'viridis',
    colorbar: bool = True,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    ax: Optional[Axes] = None,
    figsize: Tuple[float, float] = (8, 6),
) -> Axes:
    """
    Display a 2D astronomical image using astropy.visualization.

    This is a thin wrapper around astropy.visualization.imshow_norm()
    providing Heisenberg-specific defaults and convenience.

    Args:
        image: 2D numpy array to display
        wcs: Optional WCS object for coordinate transformation
        stretch: Color stretch - 'linear', 'log', 'sqrt', 'asinh'
        interval: Interval for scaling - 'minmax', 'zscale', 'percentile'
        percent: Percentile for 'percentile' interval (default 99%)
        vmin: Override minimum value for scaling
        vmax: Override maximum value for scaling
        cmap: Matplotlib colormap name
        colorbar: If True, add a colorbar
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        ax: Existing matplotlib axes (creates new figure if None)
        figsize: Figure size if creating new figure

    Returns:
        The matplotlib Axes object with the displayed image
    """
    if not HAS_ASTROPY_VIS:
        # Fallback to basic matplotlib if astropy.visualization unavailable
        return _display_image_basic(
            image, vmin=vmin, vmax=vmax, cmap=cmap,
            colorbar=colorbar, title=title, xlabel=xlabel, ylabel=ylabel,
            ax=ax, figsize=figsize
        )

    # Validate input
    if image.ndim != 2:
        raise ValueError(f"Image must be 2D, got {image.ndim}D")

    # Handle complex arrays
    if np.iscomplexobj(image):
        image = np.abs(image)

    # Create figure/axes if needed
    if ax is None:
        if wcs is not None:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111, projection=wcs)
        else:
            fig, ax = plt.subplots(figsize=figsize)

    # Build stretch
    stretch_map = {
        'linear': LinearStretch(),
        'log': LogStretch(),
        'sqrt': SqrtStretch(),
        'asinh': AsinhStretch(),
    }
    stretch_obj = stretch_map.get(stretch, LinearStretch())

    # Build interval
    if vmin is not None and vmax is not None:
        interval_obj = MinMaxInterval()
        # Override with explicit values after
    elif interval == 'zscale':
        interval_obj = ZScaleInterval()
    elif interval == 'percentile':
        interval_obj = PercentileInterval(percent)
    else:
        interval_obj = MinMaxInterval()

    # Use simple_norm for convenience
    norm = simple_norm(
        image,
        stretch=stretch,
        min_cut=vmin,
        max_cut=vmax,
    )

    # Display image
    im = ax.imshow(image, origin='lower', cmap=cmap, norm=norm)

    # Add colorbar
    if colorbar:
        plt.colorbar(im, ax=ax)

    # Format WCS axes if available
    if wcs is not None and hasattr(ax, 'coords'):
        ax.coords[0].set_axislabel(xlabel or 'Right Ascension')
        ax.coords[1].set_axislabel(ylabel or 'Declination')
    else:
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    return ax


def _display_image_basic(
    image: np.ndarray,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    cmap: str = 'viridis',
    colorbar: bool = True,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    ax: Optional[Axes] = None,
    figsize: Tuple[float, float] = (8, 6),
) -> Axes:
    """Basic matplotlib display fallback when astropy unavailable."""
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if vmin is None:
        vmin = np.nanmin(image)
    if vmax is None:
        vmax = np.nanmax(image)

    im = ax.imshow(image, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

    if colorbar:
        plt.colorbar(im, ax=ax)
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    return ax


def display_image_with_contours(
    image: np.ndarray,
    contour_levels: Optional[np.ndarray] = None,
    n_contours: int = 5,
    contour_colors: str = 'white',
    contour_linewidths: float = 0.5,
    **kwargs
) -> Axes:
    """
    Display image with contour overlays.

    Args:
        image: 2D numpy array to display
        contour_levels: Specific contour levels (if None, auto-generate)
        n_contours: Number of contours if auto-generating
        contour_colors: Color(s) for contour lines
        contour_linewidths: Line width for contours
        **kwargs: Additional arguments passed to display_image

    Returns:
        The matplotlib Axes object
    """
    ax = display_image(image, **kwargs)

    # Generate contour levels if not provided
    if contour_levels is None:
        valid_data = image[np.isfinite(image)]
        if len(valid_data) > 0:
            contour_levels = np.linspace(
                np.percentile(valid_data, 10),
                np.percentile(valid_data, 90),
                n_contours
            )
        else:
            return ax

    # Add contours
    y, x = np.mgrid[:image.shape[0], :image.shape[1]]
    ax.contour(
        x, y, image,
        levels=contour_levels,
        colors=contour_colors,
        linewidths=contour_linewidths,
    )

    return ax


def get_percentile_scaling(
    image: np.ndarray,
    low_percentile: float = 1.0,
    high_percentile: float = 99.0,
) -> Tuple[float, float]:
    """
    Get vmin/vmax from percentiles of the data.

    Args:
        image: Input image
        low_percentile: Lower percentile for vmin
        high_percentile: Upper percentile for vmax

    Returns:
        Tuple of (vmin, vmax)
    """
    valid_data = image[np.isfinite(image)]
    if len(valid_data) == 0:
        return (0.0, 1.0)

    vmin = float(np.percentile(valid_data, low_percentile))
    vmax = float(np.percentile(valid_data, high_percentile))
    return (vmin, vmax)

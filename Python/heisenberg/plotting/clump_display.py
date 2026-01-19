"""
Clump visualization functions.

Display images with clump overlays including contours, boundaries,
and peak markers.

Translates IDL clplot2d.pro functionality.
"""

from typing import List, Optional, Tuple, Union, TYPE_CHECKING

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap

from heisenberg.plotting.image_display import display_image

if TYPE_CHECKING:
    from heisenberg.peaks.clumpfind import ClumpfindResult, Clump


def plot_clumps(
    image: np.ndarray,
    clumps: 'ClumpfindResult',
    highlight_id: Optional[int] = None,
    show_boundaries: bool = True,
    show_peaks: bool = True,
    show_labels: bool = False,
    boundary_color: str = 'white',
    highlight_color: str = 'red',
    peak_color: str = 'yellow',
    peak_marker: str = '+',
    peak_size: int = 100,
    ax: Optional[Axes] = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs
) -> Axes:
    """
    Display image with clump overlays.

    Args:
        image: 2D numpy array to display
        clumps: ClumpfindResult from clumpfind2d
        highlight_id: Clump ID to highlight (1-indexed)
        show_boundaries: If True, show clump boundaries
        show_peaks: If True, mark peak positions
        show_labels: If True, add clump ID labels
        boundary_color: Color for boundary contours
        highlight_color: Color for highlighted clump
        peak_color: Color for peak markers
        peak_marker: Marker style for peaks
        peak_size: Marker size for peaks
        ax: Existing axes (creates new if None)
        figsize: Figure size if creating new figure
        **kwargs: Additional arguments passed to display_image

    Returns:
        The matplotlib Axes object
    """
    # Display base image
    ax = display_image(image, ax=ax, figsize=figsize, **kwargs)

    # Get assignment map
    assignment_map = clumps.assignment_map

    if show_boundaries:
        # Draw boundaries for each clump
        for clump in clumps.clumps:
            clump_mask = (assignment_map == clump.id).astype(float)

            # Use contour at 0.5 level to draw boundary
            color = highlight_color if clump.id == highlight_id else boundary_color
            linewidth = 2 if clump.id == highlight_id else 1

            ax.contour(
                clump_mask,
                levels=[0.5],
                colors=[color],
                linewidths=[linewidth],
            )

    if show_peaks:
        # Mark peak positions
        for clump in clumps.clumps:
            row, col = clump.peak_position
            color = highlight_color if clump.id == highlight_id else peak_color
            size = peak_size * 1.5 if clump.id == highlight_id else peak_size

            ax.scatter(
                col, row,
                marker=peak_marker,
                s=size,
                c=color,
                linewidths=2,
            )

    if show_labels:
        # Add clump ID labels
        for clump in clumps.clumps:
            row, col = clump.peak_position
            ax.annotate(
                str(clump.id),
                (col, row),
                xytext=(5, 5),
                textcoords='offset points',
                color='white',
                fontsize=8,
                fontweight='bold',
            )

    return ax


def plot_clumps_with_contours(
    image: np.ndarray,
    clumps: 'ClumpfindResult',
    contour_levels: Optional[np.ndarray] = None,
    n_contours: int = 10,
    contour_colors: str = 'gray',
    contour_alpha: float = 0.5,
    **kwargs
) -> Axes:
    """
    Display image with contour overlay and clump boundaries.

    Args:
        image: 2D numpy array to display
        clumps: ClumpfindResult from clumpfind2d
        contour_levels: Specific contour levels (auto if None)
        n_contours: Number of contours if auto
        contour_colors: Contour line color
        contour_alpha: Contour transparency
        **kwargs: Additional arguments passed to plot_clumps

    Returns:
        The matplotlib Axes object
    """
    ax = plot_clumps(image, clumps, **kwargs)

    # Generate contour levels if not provided
    if contour_levels is None:
        valid_data = image[np.isfinite(image)]
        if len(valid_data) > 0:
            contour_levels = np.linspace(
                np.percentile(valid_data, 5),
                np.percentile(valid_data, 95),
                n_contours
            )
        else:
            return ax

    # Add image contours
    ax.contour(
        image,
        levels=contour_levels,
        colors=contour_colors,
        alpha=contour_alpha,
        linewidths=0.5,
    )

    return ax


def plot_assignment_map(
    clumps: 'ClumpfindResult',
    show_colorbar: bool = True,
    cmap: str = 'tab20',
    ax: Optional[Axes] = None,
    figsize: Tuple[float, float] = (8, 6),
    title: str = 'Clump Assignment Map',
) -> Axes:
    """
    Display the clump assignment map as a colored image.

    Args:
        clumps: ClumpfindResult from clumpfind2d
        show_colorbar: If True, add colorbar
        cmap: Colormap for clump colors
        ax: Existing axes (creates new if None)
        figsize: Figure size if creating new figure
        title: Plot title

    Returns:
        The matplotlib Axes object
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Display assignment map
    im = ax.imshow(
        clumps.assignment_map,
        origin='lower',
        cmap=cmap,
        interpolation='nearest',
    )

    if show_colorbar:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Clump ID')

    ax.set_title(title)
    ax.set_xlabel('X (pixels)')
    ax.set_ylabel('Y (pixels)')

    return ax


def plot_clump_gallery(
    image: np.ndarray,
    clumps: 'ClumpfindResult',
    ncols: int = 4,
    padding: int = 5,
    figsize: Optional[Tuple[float, float]] = None,
    stretch: str = 'linear',
) -> Figure:
    """
    Display gallery of individual clump cutouts.

    Args:
        image: Original image
        clumps: ClumpfindResult from clumpfind2d
        ncols: Number of columns in gallery
        padding: Padding around clump in pixels
        figsize: Figure size (auto if None)
        stretch: Color stretch for cutouts

    Returns:
        The matplotlib Figure object
    """
    n_clumps = len(clumps.clumps)
    if n_clumps == 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, 'No clumps to display', ha='center', va='center')
        ax.set_axis_off()
        return fig

    nrows = (n_clumps + ncols - 1) // ncols

    if figsize is None:
        figsize = (3 * ncols, 3 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_2d(axes).flatten()

    for i, clump in enumerate(clumps.clumps):
        ax = axes[i]

        # Get clump bounding box
        mask = (clumps.assignment_map == clump.id)
        rows, cols = np.where(mask)

        if len(rows) == 0:
            ax.set_visible(False)
            continue

        row_min = max(0, rows.min() - padding)
        row_max = min(image.shape[0], rows.max() + padding + 1)
        col_min = max(0, cols.min() - padding)
        col_max = min(image.shape[1], cols.max() + padding + 1)

        # Extract cutout
        cutout = image[row_min:row_max, col_min:col_max]

        # Display cutout
        display_image(cutout, stretch=stretch, ax=ax, colorbar=False)

        # Mark peak
        peak_row, peak_col = clump.peak_position
        ax.scatter(
            peak_col - col_min,
            peak_row - row_min,
            marker='+',
            s=50,
            c='red',
            linewidths=1.5,
        )

        ax.set_title(f'Clump {clump.id}\n({clump.npix} pix)')

    # Hide unused axes
    for i in range(n_clumps, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    return fig


def plot_peak_comparison(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    star_peaks: List,
    gas_peaks: List,
    figsize: Tuple[float, float] = (14, 6),
    star_color: str = 'red',
    gas_color: str = 'cyan',
    stretch: str = 'log',
) -> Figure:
    """
    Side-by-side display of star and gas maps with peaks.

    Args:
        star_image: Stellar tracer image
        gas_image: Gas tracer image
        star_peaks: List of DetectedPeak for stars
        gas_peaks: List of DetectedPeak for gas
        figsize: Figure size
        star_color: Color for star peak markers
        gas_color: Color for gas peak markers
        stretch: Color stretch

    Returns:
        The matplotlib Figure object
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Star map
    display_image(star_image, stretch=stretch, ax=ax1, title='Stellar Tracer')
    for peak in star_peaks:
        ax1.scatter(peak.x, peak.y, marker='o', s=50, facecolors='none',
                    edgecolors=star_color, linewidths=1.5)

    # Gas map
    display_image(gas_image, stretch=stretch, ax=ax2, title='Gas Tracer')
    for peak in gas_peaks:
        ax2.scatter(peak.x, peak.y, marker='o', s=50, facecolors='none',
                    edgecolors=gas_color, linewidths=1.5)

    plt.tight_layout()
    return fig

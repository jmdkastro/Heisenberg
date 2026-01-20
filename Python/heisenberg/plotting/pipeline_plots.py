"""
Pipeline-integrated plotting functions.

Provides high-level plotting functions that integrate with the
Heisenberg analysis pipeline, generating all standard diagnostic
plots from analysis results.

References:
    IDL: tuningfork.pro (plotting sections)
"""

from pathlib import Path
from typing import Optional, List, Tuple, Union, TYPE_CHECKING

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

if TYPE_CHECKING:
    from heisenberg.core.tuningfork import TuningForkResult
    from heisenberg.peaks.detection import DetectedPeak


def _check_matplotlib():
    """Raise ImportError if matplotlib is not available."""
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib required for plotting: pip install matplotlib"
        )


def plot_map_with_peaks(
    image: np.ndarray,
    peaks: List['DetectedPeak'],
    ax: Optional['Axes'] = None,
    figsize: Tuple[float, float] = (8, 8),
    stretch: str = 'log',
    peak_color: str = 'red',
    peak_marker: str = 'x',
    peak_size: float = 50,
    title: Optional[str] = None,
    cmap: str = 'viridis',
    colorbar: bool = True,
    min_radius: Optional[float] = None,
    max_radius: Optional[float] = None,
    center: Optional[Tuple[float, float]] = None,
) -> Tuple['Figure', 'Axes']:
    """
    Plot an image with peak positions overlaid.

    Translates IDL plotfits + oplot peak overlay functionality.

    Args:
        image: 2D image array
        peaks: List of DetectedPeak objects
        ax: Existing axes (creates new if None)
        figsize: Figure size if creating new figure
        stretch: Color stretch ('linear', 'log', 'sqrt')
        peak_color: Color for peak markers
        peak_marker: Marker style for peaks
        peak_size: Marker size
        title: Plot title
        cmap: Colormap name
        colorbar: Whether to add colorbar
        min_radius: Inner analysis radius (optional, draws ellipse)
        max_radius: Outer analysis radius (optional, draws ellipse)
        center: Center point for radius overlays (x, y)

    Returns:
        Tuple of (figure, axes)
    """
    _check_matplotlib()
    from heisenberg.plotting.image_display import display_image

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Display image
    display_image(image, stretch=stretch, ax=ax, cmap=cmap, colorbar=colorbar)

    # Overlay peaks
    if peaks:
        x_coords = [p.x for p in peaks]
        y_coords = [p.y for p in peaks]
        ax.scatter(
            x_coords, y_coords,
            marker=peak_marker,
            s=peak_size,
            c=peak_color,
            linewidths=1.5,
            label=f'{len(peaks)} peaks',
        )

    # Draw analysis region if specified
    if min_radius is not None and center is not None:
        circle = plt.Circle(
            center, min_radius,
            fill=False, color='white', linestyle='--', linewidth=1
        )
        ax.add_patch(circle)

    if max_radius is not None and center is not None:
        circle = plt.Circle(
            center, max_radius,
            fill=False, color='white', linestyle='--', linewidth=1
        )
        ax.add_patch(circle)

    if title:
        ax.set_title(title)

    return fig, ax


def plot_star_gas_maps(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    star_peaks: List['DetectedPeak'],
    gas_peaks: List['DetectedPeak'],
    figsize: Tuple[float, float] = (14, 6),
    stretch: str = 'log',
    star_title: str = 'Stellar Tracer',
    gas_title: str = 'Gas Tracer',
    star_peak_color: str = 'red',
    gas_peak_color: str = 'cyan',
) -> 'Figure':
    """
    Side-by-side star and gas maps with peak overlays.

    Args:
        star_image: Stellar tracer image
        gas_image: Gas tracer image
        star_peaks: Detected stellar peaks
        gas_peaks: Detected gas peaks
        figsize: Figure size
        stretch: Color stretch
        star_title: Title for star panel
        gas_title: Title for gas panel
        star_peak_color: Color for stellar peak markers
        gas_peak_color: Color for gas peak markers

    Returns:
        Matplotlib Figure
    """
    _check_matplotlib()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    plot_map_with_peaks(
        star_image, star_peaks,
        ax=ax1,
        stretch=stretch,
        peak_color=star_peak_color,
        title=f'{star_title} ({len(star_peaks)} peaks)',
    )

    plot_map_with_peaks(
        gas_image, gas_peaks,
        ax=ax2,
        stretch=stretch,
        peak_color=gas_peak_color,
        title=f'{gas_title} ({len(gas_peaks)} peaks)',
    )

    plt.tight_layout()
    return fig


def plot_sensitivity_histogram(
    sensitivity: np.ndarray,
    ax: Optional['Axes'] = None,
    figsize: Tuple[float, float] = (8, 5),
    nbins: int = 50,
    title: str = 'Sensitivity Distribution',
    xlabel: str = 'Sensitivity',
    ylabel: str = 'Count',
    color: str = 'steelblue',
    log_scale: bool = False,
) -> Tuple['Figure', 'Axes']:
    """
    Plot histogram of sensitivity values.

    Translates IDL histoplot_update functionality.

    Args:
        sensitivity: 1D or 2D array of sensitivity values
        ax: Existing axes (creates new if None)
        figsize: Figure size if creating new figure
        nbins: Number of histogram bins
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        color: Histogram color
        log_scale: If True, use log scale on y-axis

    Returns:
        Tuple of (figure, axes)
    """
    _check_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Flatten and remove NaN/Inf values
    data = np.asarray(sensitivity).flatten()
    data = data[np.isfinite(data)]

    if len(data) == 0:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center',
                transform=ax.transAxes)
        ax.set_title(title)
        return fig, ax

    # Calculate bin width from data dispersion
    mean_val = np.mean(data)
    std_val = np.std(data)
    if std_val > 0:
        binwidth = std_val / (nbins / 10)
        bins = np.arange(data.min(), data.max() + binwidth, binwidth)
    else:
        bins = nbins

    ax.hist(data, bins=bins, color=color, alpha=0.7, edgecolor='black')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if log_scale:
        ax.set_yscale('log')

    # Add statistics text
    stats_text = f'Mean: {mean_val:.3g}\nStd: {std_val:.3g}\nN: {len(data)}'
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right',
            fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    return fig, ax


def plot_sensitivity_dual(
    star_sensitivity: np.ndarray,
    gas_sensitivity: np.ndarray,
    figsize: Tuple[float, float] = (12, 5),
    nbins: int = 50,
) -> 'Figure':
    """
    Side-by-side sensitivity histograms for star and gas.

    Args:
        star_sensitivity: Stellar sensitivity map
        gas_sensitivity: Gas sensitivity map
        figsize: Figure size
        nbins: Number of histogram bins

    Returns:
        Matplotlib Figure
    """
    _check_matplotlib()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    plot_sensitivity_histogram(
        star_sensitivity, ax=ax1, nbins=nbins,
        title='Stellar Sensitivity', color='salmon'
    )

    plot_sensitivity_histogram(
        gas_sensitivity, ax=ax2, nbins=nbins,
        title='Gas Sensitivity', color='lightblue'
    )

    plt.tight_layout()
    return fig


def generate_all_plots(
    result: 'TuningForkResult',
    star_image: np.ndarray,
    gas_image: np.ndarray,
    output_dir: Union[str, Path],
    galaxy_name: str = 'galaxy',
    star_sensitivity: Optional[np.ndarray] = None,
    gas_sensitivity: Optional[np.ndarray] = None,
    dpi: int = 150,
    formats: List[str] = ['png', 'pdf'],
) -> List[Path]:
    """
    Generate all standard diagnostic plots and save to files.

    This is the main integration point for pipeline plotting.
    Generates:
    - Tuning fork diagram
    - Tuning fork residuals
    - Parameter PDFs
    - Summary figure
    - Star/gas maps with peaks
    - Sensitivity histograms (if provided)

    Args:
        result: TuningForkResult from analysis
        star_image: Stellar tracer image
        gas_image: Gas tracer image
        output_dir: Directory for output files
        galaxy_name: Name for file prefixes
        star_sensitivity: Stellar sensitivity map (optional)
        gas_sensitivity: Gas sensitivity map (optional)
        dpi: Output resolution
        formats: Output file formats

    Returns:
        List of generated file paths
    """
    _check_matplotlib()
    from heisenberg.plotting.tuningfork_plot import (
        plot_tuningfork,
        plot_tuningfork_residuals,
        plot_parameter_pdfs,
        plot_summary,
        save_figure,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated_files = []

    def save_all_formats(fig, base_name):
        """Save figure in all requested formats."""
        for fmt in formats:
            path = output_dir / f'{base_name}.{fmt}'
            save_figure(fig, str(path), dpi=dpi)
            generated_files.append(path)
        plt.close(fig)

    # 1. Tuning fork diagram
    fig, _ = plot_tuningfork(result)
    save_all_formats(fig, f'{galaxy_name}_tuningfork')

    # 2. Tuning fork residuals
    fig, _ = plot_tuningfork_residuals(result)
    save_all_formats(fig, f'{galaxy_name}_residuals')

    # 3. Parameter PDFs
    fig, _ = plot_parameter_pdfs(result)
    save_all_formats(fig, f'{galaxy_name}_pdfs')

    # 4. Summary figure
    fig = plot_summary(result, suptitle=galaxy_name)
    save_all_formats(fig, f'{galaxy_name}_summary')

    # 5. Star/gas maps with peaks
    fig = plot_star_gas_maps(
        star_image, gas_image,
        result.star_peaks, result.gas_peaks,
        star_title=f'{galaxy_name} Stellar',
        gas_title=f'{galaxy_name} Gas',
    )
    save_all_formats(fig, f'{galaxy_name}_maps')

    # 6. Individual map plots
    fig, _ = plot_map_with_peaks(
        star_image, result.star_peaks,
        title=f'{galaxy_name} Stellar Tracer',
        peak_color='red',
    )
    save_all_formats(fig, f'{galaxy_name}_map_star')

    fig, _ = plot_map_with_peaks(
        gas_image, result.gas_peaks,
        title=f'{galaxy_name} Gas Tracer',
        peak_color='cyan',
    )
    save_all_formats(fig, f'{galaxy_name}_map_gas')

    # 7. Sensitivity histograms (if provided)
    if star_sensitivity is not None and gas_sensitivity is not None:
        fig = plot_sensitivity_dual(star_sensitivity, gas_sensitivity)
        save_all_formats(fig, f'{galaxy_name}_sensitivity')

    if star_sensitivity is not None:
        fig, _ = plot_sensitivity_histogram(
            star_sensitivity,
            title=f'{galaxy_name} Stellar Sensitivity',
            color='salmon',
        )
        save_all_formats(fig, f'{galaxy_name}_sensitivity_star')

    if gas_sensitivity is not None:
        fig, _ = plot_sensitivity_histogram(
            gas_sensitivity,
            title=f'{galaxy_name} Gas Sensitivity',
            color='lightblue',
        )
        save_all_formats(fig, f'{galaxy_name}_sensitivity_gas')

    return generated_files


def export_peaks_to_ds9(
    result: 'TuningForkResult',
    output_dir: Union[str, Path],
    galaxy_name: str = 'galaxy',
) -> List[Path]:
    """
    Export peak positions to DS9 region files.

    Args:
        result: TuningForkResult from analysis
        output_dir: Directory for output files
        galaxy_name: Name for file prefixes

    Returns:
        List of generated file paths
    """
    from heisenberg.plotting.ds9_regions import write_peak_regions

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated_files = []

    # Star peaks
    star_path = output_dir / f'{galaxy_name}_star_peaks.reg'
    write_peak_regions(
        result.star_peaks,
        star_path,
        color='red',
        symbol='cross',
        include_labels=True,
    )
    generated_files.append(star_path)

    # Gas peaks
    gas_path = output_dir / f'{galaxy_name}_gas_peaks.reg'
    write_peak_regions(
        result.gas_peaks,
        gas_path,
        color='cyan',
        symbol='cross',
        include_labels=True,
    )
    generated_files.append(gas_path)

    return generated_files

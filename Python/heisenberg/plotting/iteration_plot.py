"""
Iteration history plotting.

Thin wrapper around matplotlib.errorbar for tracking parameter
convergence during diffuse iteration.

Translates IDL iteration_plot.pro functionality.
"""

from typing import Optional, Tuple, Union, List

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure


def plot_iteration_history(
    iterations: np.ndarray,
    values: np.ndarray,
    errors: Optional[np.ndarray] = None,
    errors_low: Optional[np.ndarray] = None,
    errors_high: Optional[np.ndarray] = None,
    ylabel: str = 'Parameter',
    xlabel: str = 'Iteration',
    title: Optional[str] = None,
    zero_min: bool = False,
    marker: str = 'o',
    color: Optional[str] = None,
    ax: Optional[Axes] = None,
    figsize: Tuple[float, float] = (8, 5),
) -> Axes:
    """
    Plot parameter values vs iteration number with error bars.

    Supports symmetric or asymmetric error bars.

    Args:
        iterations: Array of iteration numbers (or x-values)
        values: Array of parameter values
        errors: Symmetric error bars (used if errors_low/high not provided)
        errors_low: Asymmetric lower error bounds
        errors_high: Asymmetric upper error bounds
        ylabel: Y-axis label
        xlabel: X-axis label
        title: Plot title
        zero_min: If True, set Y-axis minimum to 0
        marker: Marker style
        color: Line/marker color
        ax: Existing axes (creates new if None)
        figsize: Figure size if creating new figure

    Returns:
        The matplotlib Axes object
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Handle asymmetric errors
    if errors_low is not None and errors_high is not None:
        yerr = [errors_low, errors_high]
    elif errors is not None:
        yerr = errors
    else:
        yerr = None

    # Plot with error bars
    ax.errorbar(
        iterations,
        values,
        yerr=yerr,
        marker=marker,
        color=color,
        capsize=3,
        capthick=1,
        linewidth=1.5,
        markersize=6,
    )

    # Set labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    # Optionally set zero minimum
    if zero_min:
        ylim = ax.get_ylim()
        ax.set_ylim(0, ylim[1])

    # Integer ticks for iterations
    if np.all(iterations == iterations.astype(int)):
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    return ax


def plot_multiple_parameters(
    iterations: np.ndarray,
    parameter_dict: dict,
    errors_dict: Optional[dict] = None,
    ncols: int = 2,
    figsize: Optional[Tuple[float, float]] = None,
    zero_min: bool = False,
) -> Figure:
    """
    Plot multiple parameters vs iteration in a grid.

    Args:
        iterations: Array of iteration numbers
        parameter_dict: Dict of {param_name: values_array}
        errors_dict: Optional dict of {param_name: errors_array}
        ncols: Number of columns in grid
        figsize: Figure size (auto-calculated if None)
        zero_min: If True, set Y-axis minimum to 0 for all plots

    Returns:
        The matplotlib Figure object
    """
    n_params = len(parameter_dict)
    nrows = (n_params + ncols - 1) // ncols

    if figsize is None:
        figsize = (5 * ncols, 4 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_2d(axes).flatten()

    for i, (name, values) in enumerate(parameter_dict.items()):
        ax = axes[i]
        errors = errors_dict.get(name) if errors_dict else None

        plot_iteration_history(
            iterations,
            values,
            errors=errors,
            ylabel=name,
            title=name,
            zero_min=zero_min,
            ax=ax,
        )

    # Hide unused axes
    for i in range(n_params, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    return fig


def plot_convergence(
    iterations: np.ndarray,
    values: np.ndarray,
    tolerance: float,
    ylabel: str = 'Parameter',
    ax: Optional[Axes] = None,
    figsize: Tuple[float, float] = (8, 5),
) -> Axes:
    """
    Plot iteration history with convergence tolerance band.

    Args:
        iterations: Array of iteration numbers
        values: Array of parameter values
        tolerance: Convergence tolerance (fraction of final value)
        ylabel: Y-axis label
        ax: Existing axes (creates new if None)
        figsize: Figure size if creating new figure

    Returns:
        The matplotlib Axes object
    """
    ax = plot_iteration_history(
        iterations, values, ylabel=ylabel, ax=ax, figsize=figsize
    )

    # Add convergence band around final value
    final_value = values[-1]
    band_width = abs(final_value * tolerance)

    ax.axhline(
        final_value,
        color='green',
        linestyle='--',
        linewidth=1,
        label=f'Final: {final_value:.3g}'
    )
    ax.axhspan(
        final_value - band_width,
        final_value + band_width,
        alpha=0.2,
        color='green',
        label=f'±{tolerance*100:.1f}% tolerance'
    )

    ax.legend()
    return ax

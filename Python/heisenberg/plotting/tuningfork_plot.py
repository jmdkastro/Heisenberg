"""
Tuning fork diagram plotting for Heisenberg analysis.

Creates the characteristic "tuning fork" plot showing gas-to-stellar
flux ratios as a function of aperture size, with data and model curves.

References:
    Kruijssen & Longmore (2014), MNRAS 439, 3239
    IDL: tuningfork_plot.pro
"""

from typing import Optional, Tuple, Dict, Any
import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from heisenberg.core.tuningfork import TuningForkResult, generate_model_curve


def check_matplotlib():
    """Raise ImportError if matplotlib is not available."""
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib required for plotting: pip install matplotlib"
        )


def plot_tuningfork(
    result: TuningForkResult,
    ax: Optional["Axes"] = None,
    show_model: bool = True,
    show_data: bool = True,
    star_color: str = "C0",
    gas_color: str = "C1",
    model_n_points: int = 100,
    xlabel: str = "Aperture size (pc)",
    ylabel: str = "Gas/stellar flux ratio",
    title: Optional[str] = None,
    legend: bool = True,
    logx: bool = True,
    logy: bool = True,
) -> Tuple["Figure", "Axes"]:
    """
    Create a tuning fork diagram.

    Args:
        result: TuningForkResult from analysis
        ax: Matplotlib axes (created if None)
        show_model: Plot model curves
        show_data: Plot observed data points
        star_color: Color for stellar peak data
        gas_color: Color for gas peak data
        model_n_points: Number of points for smooth model curves
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        legend: Show legend
        logx: Logarithmic x-axis
        logy: Logarithmic y-axis

    Returns:
        Tuple of (figure, axes)
    """
    check_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure

    obs = result.observed

    # Plot observed data with error bars
    if show_data:
        ax.errorbar(
            obs.apertures,
            obs.fluxratio_star,
            yerr=obs.err_star * obs.fluxratio_star * np.log(10),
            fmt='o',
            color=star_color,
            label='Stellar peaks',
            capsize=3,
            markersize=6,
        )
        ax.errorbar(
            obs.apertures,
            obs.fluxratio_gas,
            yerr=obs.err_gas * obs.fluxratio_gas * np.log(10),
            fmt='s',
            color=gas_color,
            label='Gas peaks',
            capsize=3,
            markersize=6,
        )

    # Plot model curves
    if show_model:
        apertures, model_star, model_gas = generate_model_curve(
            result.fit, n_points=model_n_points
        )
        ax.plot(
            apertures,
            model_star,
            '-',
            color=star_color,
            alpha=0.8,
            linewidth=2,
            label='Model (stellar)' if show_data else 'Stellar peaks',
        )
        ax.plot(
            apertures,
            model_gas,
            '-',
            color=gas_color,
            alpha=0.8,
            linewidth=2,
            label='Model (gas)' if show_data else 'Gas peaks',
        )

    # Horizontal line at ratio = 1
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5, linewidth=1)

    # Axis settings
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    if legend:
        ax.legend(loc='best')

    ax.grid(True, alpha=0.3)

    return fig, ax


def plot_tuningfork_residuals(
    result: TuningForkResult,
    ax: Optional["Axes"] = None,
    star_color: str = "C0",
    gas_color: str = "C1",
    xlabel: str = "Aperture size (pc)",
    ylabel: str = "Residual (data - model)",
) -> Tuple["Figure", "Axes"]:
    """
    Plot residuals between observed data and model.

    Args:
        result: TuningForkResult from analysis
        ax: Matplotlib axes (created if None)
        star_color: Color for stellar peak residuals
        gas_color: Color for gas peak residuals
        xlabel: X-axis label
        ylabel: Y-axis label

    Returns:
        Tuple of (figure, axes)
    """
    check_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.figure

    obs = result.observed

    # Get model values at observed apertures
    apertures, model_star, model_gas = generate_model_curve(
        result.fit, apertures=obs.apertures
    )

    # Compute residuals
    resid_star = obs.fluxratio_star - model_star
    resid_gas = obs.fluxratio_gas - model_gas

    # Plot residuals
    ax.errorbar(
        obs.apertures,
        resid_star,
        yerr=obs.err_star * obs.fluxratio_star * np.log(10),
        fmt='o',
        color=star_color,
        label='Stellar peaks',
        capsize=3,
    )
    ax.errorbar(
        obs.apertures,
        resid_gas,
        yerr=obs.err_gas * obs.fluxratio_gas * np.log(10),
        fmt='s',
        color=gas_color,
        label='Gas peaks',
        capsize=3,
    )

    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    return fig, ax


def plot_parameter_pdfs(
    result: TuningForkResult,
    figsize: Tuple[float, float] = (12, 4),
) -> Tuple["Figure", Tuple["Axes", "Axes", "Axes"]]:
    """
    Plot probability density functions for fitted parameters.

    Args:
        result: TuningForkResult from analysis
        figsize: Figure size

    Returns:
        Tuple of (figure, (ax_tgas, ax_tover, ax_lambda))
    """
    check_matplotlib()

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    ax_tgas, ax_tover, ax_lambda = axes

    fit = result.fit

    # t_gas PDF
    ax_tgas.plot(fit.tgas_arr, fit.prob_tgas, 'b-', linewidth=2)
    ax_tgas.axvline(fit.tgas, color='r', linestyle='--', label='Best fit')
    ax_tgas.fill_between(fit.tgas_arr, fit.prob_tgas, alpha=0.3)
    ax_tgas.set_xlabel(r'$t_{\rm gas}$ (Myr)')
    ax_tgas.set_ylabel('Probability density')
    ax_tgas.set_title(
        f'$t_{{\\rm gas}} = {fit.tgas:.1f}^{{+{fit.tgas_errmax:.1f}}}_{{-{fit.tgas_errmin:.1f}}}$ Myr'
    )

    # t_over PDF
    ax_tover.plot(fit.tover_arr, fit.prob_tover, 'g-', linewidth=2)
    ax_tover.axvline(fit.tover, color='r', linestyle='--', label='Best fit')
    ax_tover.fill_between(fit.tover_arr, fit.prob_tover, alpha=0.3, color='g')
    ax_tover.set_xlabel(r'$t_{\rm over}$ (Myr)')
    ax_tover.set_ylabel('Probability density')
    ax_tover.set_title(
        f'$t_{{\\rm over}} = {fit.tover:.2f}^{{+{fit.tover_errmax:.2f}}}_{{-{fit.tover_errmin:.2f}}}$ Myr'
    )

    # lambda PDF
    ax_lambda.plot(fit.lambda_arr, fit.prob_lambda, 'm-', linewidth=2)
    ax_lambda.axvline(fit.lambda_, color='r', linestyle='--', label='Best fit')
    ax_lambda.fill_between(fit.lambda_arr, fit.prob_lambda, alpha=0.3, color='m')
    ax_lambda.set_xlabel(r'$\lambda$ (pc)')
    ax_lambda.set_ylabel('Probability density')
    ax_lambda.set_title(
        f'$\\lambda = {fit.lambda_:.0f}^{{+{fit.lambda_errmax:.0f}}}_{{-{fit.lambda_errmin:.0f}}}$ pc'
    )

    plt.tight_layout()
    return fig, (ax_tgas, ax_tover, ax_lambda)


def plot_summary(
    result: TuningForkResult,
    figsize: Tuple[float, float] = (14, 10),
    suptitle: Optional[str] = None,
) -> "Figure":
    """
    Create a multi-panel summary figure.

    Includes tuning fork diagram, residuals, and parameter PDFs.

    Args:
        result: TuningForkResult from analysis
        figsize: Figure size
        suptitle: Overall figure title

    Returns:
        Matplotlib figure
    """
    check_matplotlib()

    fig = plt.figure(figsize=figsize)

    # Create grid
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # Tuning fork diagram (top, spans 2 columns)
    ax_tf = fig.add_subplot(gs[0:2, 0:2])
    plot_tuningfork(result, ax=ax_tf)

    # Residuals (middle right)
    ax_resid = fig.add_subplot(gs[0, 2])
    plot_tuningfork_residuals(result, ax=ax_resid)
    ax_resid.set_title('Residuals')

    # PDFs (bottom row)
    fit = result.fit

    ax_tgas = fig.add_subplot(gs[2, 0])
    ax_tgas.plot(fit.tgas_arr, fit.prob_tgas, 'b-', linewidth=2)
    ax_tgas.axvline(fit.tgas, color='r', linestyle='--')
    ax_tgas.fill_between(fit.tgas_arr, fit.prob_tgas, alpha=0.3)
    ax_tgas.set_xlabel(r'$t_{\rm gas}$ (Myr)')
    ax_tgas.set_title(f'$t_{{\\rm gas}} = {fit.tgas:.1f}$ Myr')

    ax_tover = fig.add_subplot(gs[2, 1])
    ax_tover.plot(fit.tover_arr, fit.prob_tover, 'g-', linewidth=2)
    ax_tover.axvline(fit.tover, color='r', linestyle='--')
    ax_tover.fill_between(fit.tover_arr, fit.prob_tover, alpha=0.3, color='g')
    ax_tover.set_xlabel(r'$t_{\rm over}$ (Myr)')
    ax_tover.set_title(f'$t_{{\\rm over}} = {fit.tover:.2f}$ Myr')

    ax_lambda = fig.add_subplot(gs[2, 2])
    ax_lambda.plot(fit.lambda_arr, fit.prob_lambda, 'm-', linewidth=2)
    ax_lambda.axvline(fit.lambda_, color='r', linestyle='--')
    ax_lambda.fill_between(fit.lambda_arr, fit.prob_lambda, alpha=0.3, color='m')
    ax_lambda.set_xlabel(r'$\lambda$ (pc)')
    ax_lambda.set_title(f'$\\lambda = {fit.lambda_:.0f}$ pc')

    # Text summary (middle right bottom)
    ax_text = fig.add_subplot(gs[1, 2])
    ax_text.axis('off')

    derived = result.derived
    text_lines = [
        'Derived quantities:',
        f"$t_{{\\rm total}} = {derived.get('ttotal', 0):.1f}$ Myr",
        f"$\\epsilon_{{\\rm sf}} = {derived.get('esf', 0):.1%}$",
        f"$v_{{\\rm fb}} = {derived.get('vfb', 0):.1f}$ km/s",
        '',
        f"$\\chi^2_{{\\rm min}} = {fit.chi2_min:.2f}$",
        f"ndof = {fit.ndof}",
    ]
    ax_text.text(
        0.1, 0.9, '\n'.join(text_lines),
        transform=ax_text.transAxes,
        fontsize=10,
        verticalalignment='top',
        fontfamily='monospace',
    )

    if suptitle:
        fig.suptitle(suptitle, fontsize=14, fontweight='bold')

    return fig


def save_figure(
    fig: "Figure",
    path: str,
    dpi: int = 150,
    bbox_inches: str = 'tight',
    **kwargs,
) -> None:
    """
    Save figure to file.

    Args:
        fig: Matplotlib figure
        path: Output file path
        dpi: Resolution
        bbox_inches: Bounding box setting
        **kwargs: Additional arguments to savefig
    """
    check_matplotlib()
    fig.savefig(path, dpi=dpi, bbox_inches=bbox_inches, **kwargs)

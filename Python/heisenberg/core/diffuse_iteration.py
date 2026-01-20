"""
Diffuse Iteration Controller.

Implements iterative Fourier filtering for removing diffuse emission,
running the Heisenberg analysis repeatedly until convergence.

The iteration loop:
1. Run tuningfork analysis to get lambda
2. Apply Fourier high-pass filter at lambda scale
3. Re-run analysis on filtered images
4. Check convergence criterion
5. Repeat until converged or max iterations reached

References:
    IDL: diffuse_iteration.pro
"""

import shutil
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass, field

import numpy as np

from heisenberg.config import HeisenbergConfig
from heisenberg.imaging import read_fits, write_fits, get_platescale
from heisenberg.fourier import apply_fourier_filter
from heisenberg.core.tuningfork import (
    run_tuningfork,
    TuningForkResult,
    TuningForkConfig,
)
from heisenberg.plotting import plot_iteration_history


@dataclass
class IterationHistory:
    """
    Track parameter values across iterations.

    Attributes:
        tgas: Gas phase timescale history (Myr)
        tgas_errmin: Lower error on tgas
        tgas_errmax: Upper error on tgas
        tover: Overlap timescale history (Myr)
        tover_errmin: Lower error on tover
        tover_errmax: Upper error on tover
        lambda_: Mean separation history (pc)
        lambda_errmin: Lower error on lambda
        lambda_errmax: Upper error on lambda
        npixmin: Peak detection npixmin history
        nsigma: Peak detection nsigma history
        converged: Whether iteration converged
        n_iterations: Number of iterations completed
    """
    tgas: List[float] = field(default_factory=list)
    tgas_errmin: List[float] = field(default_factory=list)
    tgas_errmax: List[float] = field(default_factory=list)
    tover: List[float] = field(default_factory=list)
    tover_errmin: List[float] = field(default_factory=list)
    tover_errmax: List[float] = field(default_factory=list)
    lambda_: List[float] = field(default_factory=list)
    lambda_errmin: List[float] = field(default_factory=list)
    lambda_errmax: List[float] = field(default_factory=list)
    npixmin: List[int] = field(default_factory=list)
    nsigma: List[float] = field(default_factory=list)
    converged: bool = False
    n_iterations: int = 0

    def append_result(
        self,
        result: TuningForkResult,
        npixmin: int,
        nsigma: float,
    ) -> None:
        """Append a tuningfork result to the history."""
        fit = result.fit
        self.tgas.append(fit.tgas)
        self.tgas_errmin.append(fit.tgas_errmin)
        self.tgas_errmax.append(fit.tgas_errmax)
        self.tover.append(fit.tover)
        self.tover_errmin.append(fit.tover_errmin)
        self.tover_errmax.append(fit.tover_errmax)
        self.lambda_.append(fit.lambda_)
        self.lambda_errmin.append(fit.lambda_errmin)
        self.lambda_errmax.append(fit.lambda_errmax)
        self.npixmin.append(npixmin)
        self.nsigma.append(nsigma)
        self.n_iterations = len(self.tgas)


@dataclass
class DiffuseIterationConfig:
    """
    Configuration for diffuse iteration.

    Attributes:
        use_guess: If True, apply initial filter at initial_guess scale
        initial_guess: Initial lambda guess for first filter (pc)
        iter_criterion: Convergence criterion (fractional change in lambda)
        iter_crit_len: Number of consecutive iterations that must satisfy criterion
        iter_nmax: Maximum number of iterations
        iter_filter: Filter type (0=butterworth, 1=gaussian, 2=ideal)
        iter_bwo: Butterworth order (if using butterworth filter)
        iter_len_conv: Filter length conversion factor (filter at lambda * iter_len_conv)
        use_noisecut: If True, apply noise threshold after filtering
        noisethresh_s: Noise threshold for stellar image
        noisethresh_g: Noise threshold for gas image
        zero_negatives: If True, set negative values to zero after filtering
    """
    use_guess: bool = False
    initial_guess: float = 200.0
    iter_criterion: float = 0.05
    iter_crit_len: int = 2
    iter_nmax: int = 10
    iter_filter: int = 2  # 0=butterworth, 1=gaussian, 2=ideal
    iter_bwo: int = 2
    iter_len_conv: float = 2.0
    use_noisecut: bool = False
    noisethresh_s: float = 0.0
    noisethresh_g: float = 0.0
    zero_negatives: bool = True


def apply_image_threshold(
    image: np.ndarray,
    threshold: Optional[float] = None,
    zero_negatives: bool = True,
) -> np.ndarray:
    """
    Apply threshold to image, zeroing values below threshold.

    This is a Python translation of IDL apply_image_threshold.pro.

    Args:
        image: Image array (modified in place)
        threshold: Values below this are set to zero
        zero_negatives: If True, also zero all negative values

    Returns:
        Modified image array
    """
    result = image.copy()

    if threshold is not None and threshold >= 0:
        result[result < threshold] = 0.0

    if zero_negatives:
        result[result < 0] = 0.0

    return result


def iter_image_postprocess(
    image: np.ndarray,
    threshold: Optional[float] = None,
    zero_negatives: bool = True,
) -> np.ndarray:
    """
    Post-process filtered image for iteration.

    This is a Python translation of IDL iter_image_postprocess.pro.

    Args:
        image: Image array to process
        threshold: Noise threshold
        zero_negatives: If True, zero negative values

    Returns:
        Processed image array
    """
    return apply_image_threshold(image, threshold, zero_negatives)


def _check_convergence(
    history: IterationHistory,
    criterion: float,
    crit_len: int,
) -> bool:
    """
    Check if lambda has converged.

    Convergence is satisfied when the fractional change in lambda
    is less than the criterion for crit_len consecutive iterations.

    Args:
        history: Iteration history
        criterion: Convergence criterion (fractional)
        crit_len: Number of consecutive iterations to check

    Returns:
        True if converged
    """
    if len(history.lambda_) <= crit_len:
        return False

    lambda_current = history.lambda_[-1]

    for i in range(1, crit_len + 1):
        lambda_prev = history.lambda_[-(i + 1)]
        frac_change = abs((lambda_current - lambda_prev) / lambda_current)
        if frac_change >= criterion:
            return False

    return True


def _get_filter_type(filter_int: int) -> str:
    """Convert integer filter type to string."""
    filter_map = {0: 'butterworth', 1: 'gaussian', 2: 'ideal'}
    return filter_map.get(filter_int, 'ideal')


def run_diffuse_iteration(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    star_pix_to_pc: float,
    gas_pix_to_pc: float,
    config: Optional[DiffuseIterationConfig] = None,
    tuningfork_config: Optional[TuningForkConfig] = None,
    star_sensitivity: Optional[np.ndarray] = None,
    gas_sensitivity: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
    output_dir: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[TuningForkResult, IterationHistory]:
    """
    Run iterative diffuse filtering analysis.

    This is the main entry point for iterative Heisenberg analysis,
    implementing the algorithm from IDL diffuse_iteration.pro.

    The process:
    1. Optionally apply initial guess filter
    2. Run tuningfork analysis to get lambda
    3. Apply high-pass Fourier filter at lambda scale
    4. Re-run analysis on filtered images
    5. Check convergence (fractional change in lambda < criterion)
    6. Repeat until converged or max iterations reached

    Args:
        star_image: Original stellar tracer image
        gas_image: Original gas tracer image
        star_pix_to_pc: Pixel scale for stellar image (pc/pixel)
        gas_pix_to_pc: Pixel scale for gas image (pc/pixel)
        config: DiffuseIterationConfig with iteration parameters
        tuningfork_config: TuningForkConfig for each analysis
        star_sensitivity: Sensitivity map for stellar image
        gas_sensitivity: Sensitivity map for gas image
        mask: Optional mask (1=valid, 0=masked)
        output_dir: Directory for iteration outputs (optional)
        verbose: If True, print progress messages

    Returns:
        Tuple of (final TuningForkResult, IterationHistory)
    """
    if config is None:
        config = DiffuseIterationConfig()

    if tuningfork_config is None:
        tuningfork_config = TuningForkConfig()

    # Set up output directory
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        iter_plotdir = output_dir / 'iteration_plots'
        iter_plotdir.mkdir(exist_ok=True)
        iter_datadir = output_dir / 'iteration_data'
        iter_datadir.mkdir(exist_ok=True)

    # Store base images
    base_star = star_image.copy()
    base_gas = gas_image.copy()

    # Current working images
    current_star = star_image.copy()
    current_gas = gas_image.copy()

    history = IterationHistory()
    filter_type = _get_filter_type(config.iter_filter)

    # Apply initial guess filtering if requested
    if config.use_guess:
        if verbose:
            print(f"=Iter==> Applying initial guess filter at {config.initial_guess:.1f} pc")

        # Filter star image
        star_cut_pix = (config.initial_guess / star_pix_to_pc) * config.iter_len_conv
        current_star = apply_fourier_filter(
            base_star,
            cutoff=star_cut_pix,
            filter_type=filter_type,
            pass_type='high',
            order=config.iter_bwo,
        )
        current_star = iter_image_postprocess(
            current_star,
            threshold=config.noisethresh_s if config.use_noisecut else None,
            zero_negatives=config.zero_negatives,
        )

        # Filter gas image
        gas_cut_pix = (config.initial_guess / gas_pix_to_pc) * config.iter_len_conv
        current_gas = apply_fourier_filter(
            base_gas,
            cutoff=gas_cut_pix,
            filter_type=filter_type,
            pass_type='high',
            order=config.iter_bwo,
        )
        current_gas = iter_image_postprocess(
            current_gas,
            threshold=config.noisethresh_g if config.use_noisecut else None,
            zero_negatives=config.zero_negatives,
        )

    # Iteration loop
    result = None
    for iter_num in range(config.iter_nmax):
        if verbose:
            print(f"=Iter==> Running iteration {iter_num}")

        # Run tuningfork analysis
        try:
            result = run_tuningfork(
                star_image=current_star,
                gas_image=current_gas,
                pixel_scale=star_pix_to_pc,  # Use star pixel scale
                config=tuningfork_config,
                star_sensitivity=star_sensitivity,
                gas_sensitivity=gas_sensitivity,
                mask=mask,
            )
        except ValueError as e:
            if verbose:
                print(f"=Iter==> Analysis failed: {e}")
            break

        # Record results
        history.append_result(
            result,
            npixmin=tuningfork_config.npixmin_star,
            nsigma=tuningfork_config.nsigma_star,
        )

        lambda_val = result.fit.lambda_
        if verbose:
            print(f"=Iter==> Iteration {iter_num}: lambda = {lambda_val:.1f} pc, "
                  f"tgas = {result.fit.tgas:.2f} Myr, tover = {result.fit.tover:.2f} Myr")

        # Save iteration output if output_dir specified
        if output_dir is not None:
            _save_iteration_output(output_dir, iter_num, result, current_star, current_gas)

        # Check convergence
        if _check_convergence(history, config.iter_criterion, config.iter_crit_len):
            if verbose:
                print("=Iter==> Convergence criterion satisfied!")
            history.converged = True
            break

        # Apply Fourier filter for next iteration
        if iter_num < config.iter_nmax - 1:
            # Filter star image
            star_cut_pix = (lambda_val / star_pix_to_pc) * config.iter_len_conv
            current_star = apply_fourier_filter(
                base_star,
                cutoff=star_cut_pix,
                filter_type=filter_type,
                pass_type='high',
                order=config.iter_bwo,
            )
            current_star = iter_image_postprocess(
                current_star,
                threshold=config.noisethresh_s if config.use_noisecut else None,
                zero_negatives=config.zero_negatives,
            )

            # Filter gas image
            gas_cut_pix = (lambda_val / gas_pix_to_pc) * config.iter_len_conv
            current_gas = apply_fourier_filter(
                base_gas,
                cutoff=gas_cut_pix,
                filter_type=filter_type,
                pass_type='high',
                order=config.iter_bwo,
            )
            current_gas = iter_image_postprocess(
                current_gas,
                threshold=config.noisethresh_g if config.use_noisecut else None,
                zero_negatives=config.zero_negatives,
            )

    # Generate iteration plots if output_dir specified
    if output_dir is not None and len(history.tgas) > 1:
        _generate_iteration_plots(output_dir / 'iteration_plots', history)

    if result is None:
        raise RuntimeError("Diffuse iteration failed to complete any successful runs")

    return result, history


def _save_iteration_output(
    output_dir: Path,
    iter_num: int,
    result: TuningForkResult,
    star_image: np.ndarray,
    gas_image: np.ndarray,
) -> None:
    """Save outputs for a single iteration."""
    iter_datadir = output_dir / 'iteration_data'

    # Save filtered images
    # Note: Would need header info for proper FITS writing
    # For now, save as numpy arrays
    np.save(iter_datadir / f'starfile_iter{iter_num}.npy', star_image)
    np.save(iter_datadir / f'gasfile_iter{iter_num}.npy', gas_image)

    # Save fit results
    fit = result.fit
    with open(iter_datadir / f'result_iter{iter_num}.txt', 'w') as f:
        f.write(f"# Iteration {iter_num} results\n")
        f.write(f"tgas {fit.tgas:.6f}\n")
        f.write(f"tgas_errmin {fit.tgas_errmin:.6f}\n")
        f.write(f"tgas_errmax {fit.tgas_errmax:.6f}\n")
        f.write(f"tover {fit.tover:.6f}\n")
        f.write(f"tover_errmin {fit.tover_errmin:.6f}\n")
        f.write(f"tover_errmax {fit.tover_errmax:.6f}\n")
        f.write(f"lambda {fit.lambda_:.6f}\n")
        f.write(f"lambda_errmin {fit.lambda_errmin:.6f}\n")
        f.write(f"lambda_errmax {fit.lambda_errmax:.6f}\n")


def _generate_iteration_plots(
    plot_dir: Path,
    history: IterationHistory,
) -> None:
    """Generate iteration tracking plots."""
    import matplotlib.pyplot as plt

    iterations = np.arange(history.n_iterations)

    # Plot each tracked parameter
    params = [
        ('tgas', 't_gas [Myr]', history.tgas, history.tgas_errmin, history.tgas_errmax),
        ('tover', 't_over [Myr]', history.tover, history.tover_errmin, history.tover_errmax),
        ('lambda', 'λ [pc]', history.lambda_, history.lambda_errmin, history.lambda_errmax),
    ]

    for name, label, values, err_low, err_high in params:
        fig, ax = plt.subplots(figsize=(8, 5))

        values_arr = np.array(values)
        err_low_arr = np.array(err_low)
        err_high_arr = np.array(err_high)

        plot_iteration_history(
            iterations,
            values_arr,
            errors_low=err_low_arr,
            errors_high=err_high_arr,
            ylabel=label,
            xlabel='Iteration',
            title=f'{label} vs Iteration',
            ax=ax,
        )

        fig.tight_layout()
        fig.savefig(plot_dir / f'{name}_iteration.png', dpi=150)
        plt.close(fig)

    # Simple parameter plots (no errors)
    simple_params = [
        ('npixmin', 'N_pixmin', history.npixmin),
        ('nsigma', 'N_sigma', history.nsigma),
    ]

    for name, label, values in simple_params:
        fig, ax = plt.subplots(figsize=(8, 5))

        plot_iteration_history(
            iterations,
            np.array(values),
            ylabel=label,
            xlabel='Iteration',
            title=f'{label} vs Iteration',
            ax=ax,
        )

        fig.tight_layout()
        fig.savefig(plot_dir / f'{name}_iteration.png', dpi=150)
        plt.close(fig)


def write_iteration_report(
    filepath: Path,
    history: IterationHistory,
) -> None:
    """
    Write iteration history to a report file.

    Args:
        filepath: Output file path
        history: IterationHistory to write
    """
    with open(filepath, 'w') as f:
        # Header
        f.write("# Heisenberg Diffuse Iteration Report\n")
        f.write(f"# Converged: {history.converged}\n")
        f.write(f"# N_iterations: {history.n_iterations}\n")
        f.write("#\n")
        f.write("# iter  tgas[Myr]  tgas_emin  tgas_emax  "
                "tover[Myr]  tover_emin  tover_emax  "
                "lambda[pc]  lambda_emin  lambda_emax\n")

        for i in range(history.n_iterations):
            f.write(f"{i:4d}  "
                    f"{history.tgas[i]:10.4f}  {history.tgas_errmin[i]:10.4f}  "
                    f"{history.tgas_errmax[i]:10.4f}  "
                    f"{history.tover[i]:10.4f}  {history.tover_errmin[i]:10.4f}  "
                    f"{history.tover_errmax[i]:10.4f}  "
                    f"{history.lambda_[i]:10.2f}  {history.lambda_errmin[i]:10.2f}  "
                    f"{history.lambda_errmax[i]:10.2f}\n")

"""
Tuning Fork Analysis Pipeline.

Main orchestration module for the Heisenberg uncertainty principle analysis.
Coordinates peak identification, aperture photometry, flux ratio calculation,
and KL14 uncertainty principle fitting.

The "tuning fork" diagram shows gas-to-stellar flux ratios as a function
of aperture size, with characteristic divergence between stellar-peak and
gas-peak measurements at small scales.

References:
    Kruijssen & Longmore (2014), MNRAS 439, 3239
    Kruijssen et al. (2018), MNRAS 479, 1866
    IDL: tuningfork.pro
"""

import numpy as np
from typing import Optional, List, Tuple, Dict, Any, Union
from dataclasses import dataclass, field
from pathlib import Path

from heisenberg.imaging import (
    generate_aperture_sizes,
    measure_aperture_flux,
    ApertureFluxResult,
)
from heisenberg.imaging.aperture import monte_carlo_flux_ratios
from heisenberg.peaks import (
    find_peaks,
    PeakDetectionConfig,
    DetectedPeak,
    peaks_to_coords,
)
from heisenberg.core.fitting import fit_tuningfork, FitResult
from heisenberg.core.derived import (
    f_esf,
    f_vfb,
    f_etainst,
    f_etaavg,
)
from heisenberg.peaks.nearest_neighbour import (
    nearest_neighbour_from_peaks,
    NearestNeighbourResult,
)


@dataclass
class TuningForkData:
    """
    Observed data for tuning fork diagram.

    Attributes:
        apertures: Array of aperture sizes (pc)
        fluxratio_star: Gas/stellar flux ratio at stellar peaks
        fluxratio_gas: Gas/stellar flux ratio at gas peaks
        err_star: Errors on stellar flux ratios
        err_gas: Errors on gas flux ratios
        n_star_peaks: Number of stellar peaks used
        n_gas_peaks: Number of gas peaks used
    """
    apertures: np.ndarray
    fluxratio_star: np.ndarray
    fluxratio_gas: np.ndarray
    err_star: np.ndarray
    err_gas: np.ndarray
    n_star_peaks: np.ndarray
    n_gas_peaks: np.ndarray


@dataclass
class TuningForkResult:
    """
    Complete result from tuning fork analysis.

    Attributes:
        observed: Observed tuning fork data
        fit: FitResult from KL14 uncertainty principle fitting
        star_peaks: List of detected stellar peaks
        gas_peaks: List of detected gas peaks
        derived: Dictionary of derived physical quantities
        config: Configuration used for analysis
        star_nn: Nearest-neighbour statistics for stellar peaks
        gas_nn: Nearest-neighbour statistics for gas peaks
    """
    observed: TuningForkData
    fit: FitResult
    star_peaks: List[DetectedPeak]
    gas_peaks: List[DetectedPeak]
    derived: Dict[str, float] = field(default_factory=dict)
    config: Dict[str, Any] = field(default_factory=dict)
    star_nn: Optional[NearestNeighbourResult] = None
    gas_nn: Optional[NearestNeighbourResult] = None


@dataclass
class TuningForkConfig:
    """
    Configuration for tuning fork analysis.

    Attributes:
        lap_min: Minimum aperture size (pc)
        lap_max: Maximum aperture size (pc)
        n_apertures: Number of aperture sizes
        n_mc: Number of Monte Carlo realizations
        tstar: Stellar tracer visibility time (Myr)
        peak_prof: Peak profile model (0=point, 1=disc, 2=Gaussian)
        npixmin_star: Minimum pixels per stellar peak
        npixmin_gas: Minimum pixels per gas peak
        nsigma_star: Detection threshold for stellar peaks
        nsigma_gas: Detection threshold for gas peaks
        min_area_frac: Minimum fractional area for valid measurement
        seed: Random seed for reproducibility
        tstar_incl: If True, tstar includes overlap phase
        tgas_min: Minimum tgas for fitting (Myr)
        tgas_max: Maximum tgas for fitting (Myr)
        tover_min: Minimum tover for fitting (Myr)
        ndepth: Number of fitting refinement iterations
        ntry: Grid points per dimension in fitting
    """
    lap_min: float = 50.0
    lap_max: float = 1000.0
    n_apertures: int = 15
    n_mc: int = 100
    tstar: float = 10.0
    peak_prof: int = 2
    npixmin_star: int = 20
    npixmin_gas: int = 20
    nsigma_star: float = 5.0
    nsigma_gas: float = 5.0
    min_area_frac: float = 0.5
    seed: Optional[int] = None
    tstar_incl: bool = False
    tgas_min: float = 0.5
    tgas_max: float = 50.0
    tover_min: float = 0.1
    ndepth: int = 4
    ntry: int = 51


def run_tuningfork(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    pixel_scale: float,
    config: Optional[TuningForkConfig] = None,
    star_sensitivity: Optional[np.ndarray] = None,
    gas_sensitivity: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
    star_beam_fwhm: float = 0.0,
    gas_beam_fwhm: float = 0.0,
) -> TuningForkResult:
    """
    Run the complete tuning fork analysis pipeline.

    This is the main entry point for the Heisenberg analysis. It:
    1. Identifies peaks in stellar and gas maps
    2. Measures flux at peaks across multiple aperture sizes
    3. Computes flux ratios with Monte Carlo uncertainty
    4. Fits the KL14 uncertainty principle to derive timescales

    Args:
        star_image: Stellar tracer map (e.g., H-alpha, FUV)
        gas_image: Gas tracer map (e.g., CO)
        pixel_scale: Physical size per pixel (e.g., pc/pixel)
        config: TuningForkConfig with analysis parameters
        star_sensitivity: Sensitivity map for stellar image
        gas_sensitivity: Sensitivity map for gas image
        mask: Optional mask (1=valid, 0=masked)
        star_beam_fwhm: Current beam FWHM for stellar image (in pixels)
        gas_beam_fwhm: Current beam FWHM for gas image (in pixels)

    Returns:
        TuningForkResult with observed data, fit results, and derived quantities
    """
    if config is None:
        config = TuningForkConfig()

    # Step 1: Generate aperture sizes
    apertures = generate_aperture_sizes(
        config.lap_min,
        config.lap_max,
        config.n_apertures,
    )

    # Step 2: Detect peaks
    star_peak_config = PeakDetectionConfig(
        npixmin=config.npixmin_star,
        nsigma=config.nsigma_star,
    )
    gas_peak_config = PeakDetectionConfig(
        npixmin=config.npixmin_gas,
        nsigma=config.nsigma_gas,
    )

    star_peaks = find_peaks(
        star_image,
        config=star_peak_config,
        sensitivity=star_sensitivity,
    )
    gas_peaks = find_peaks(
        gas_image,
        config=gas_peak_config,
        sensitivity=gas_sensitivity,
    )

    if len(star_peaks) < 2:
        raise ValueError(f"Insufficient stellar peaks: {len(star_peaks)} (need >= 2)")
    if len(gas_peaks) < 2:
        raise ValueError(f"Insufficient gas peaks: {len(gas_peaks)} (need >= 2)")

    # Step 2b: Calculate nearest-neighbour statistics (for overlap correction)
    star_nn = nearest_neighbour_from_peaks(star_peaks, mask=mask, seed=config.seed)
    gas_nn = nearest_neighbour_from_peaks(gas_peaks, mask=mask, seed=config.seed)

    # Step 3: Measure aperture flux
    star_coords = peaks_to_coords(star_peaks)
    gas_coords = peaks_to_coords(gas_peaks)

    star_flux_result, gas_flux_result = measure_aperture_flux(
        star_image=star_image,
        gas_image=gas_image,
        star_peaks=star_coords,
        gas_peaks=gas_coords,
        apertures=apertures / pixel_scale,  # Convert to pixels
        pixel_scale=pixel_scale,
        star_beam_fwhm=star_beam_fwhm,
        gas_beam_fwhm=gas_beam_fwhm,
        mask=mask,
    )

    # Step 4: Compute flux ratios with Monte Carlo
    mean_star, std_star, n_star = monte_carlo_flux_ratios(
        star_flux_result,
        n_mc=config.n_mc,
        min_area_frac=config.min_area_frac,
        seed=config.seed,
    )

    mean_gas, std_gas, n_gas = monte_carlo_flux_ratios(
        gas_flux_result,
        n_mc=config.n_mc,
        min_area_frac=config.min_area_frac,
        seed=config.seed,
    )

    # Convert to log errors
    err_star_log = np.where(mean_star > 0, std_star / mean_star / np.log(10), 0.1)
    err_gas_log = np.where(mean_gas > 0, std_gas / mean_gas / np.log(10), 0.1)

    # Create observed data structure
    observed = TuningForkData(
        apertures=apertures,
        fluxratio_star=mean_star,
        fluxratio_gas=mean_gas,
        err_star=err_star_log,
        err_gas=err_gas_log,
        n_star_peaks=n_star,
        n_gas_peaks=n_gas,
    )

    # Step 5: Fit KL14 uncertainty principle
    # Create default beta and overlap arrays (unity for simple case)
    # beta_star/beta_gas are interpolation tables: value at each fstarover/fgasover
    n_ap = len(apertures)
    n_interp = 100
    fstarover = np.linspace(0, 1, n_interp)
    fgasover = np.linspace(0, 1, n_interp)
    beta_star = np.ones(n_interp)  # Same length as fstarover
    beta_gas = np.ones(n_interp)   # Same length as fgasover
    surfcontrasts = np.ones(n_ap) * 5.0  # Default surface contrast
    surfcontrastg = np.ones(n_ap) * 5.0

    fit = fit_tuningfork(
        fluxratio_star=mean_star,
        fluxratio_gas=mean_gas,
        err_star_log=err_star_log,
        err_gas_log=err_gas_log,
        tstariso=config.tstar,
        beta_star=beta_star,
        beta_gas=beta_gas,
        fstarover=fstarover,
        fgasover=fgasover,
        apertures_star=apertures,
        apertures_gas=apertures,
        surfcontrasts=surfcontrasts,
        surfcontrastg=surfcontrastg,
        peak_prof=config.peak_prof,
        tstar_incl=config.tstar_incl,
        tgasmini=config.tgas_min,
        tgasmaxi=config.tgas_max,
        tovermini=config.tover_min,
        ndepth=config.ndepth,
        ntry=config.ntry,
        nfitstar=n_star,
        nfitgas=n_gas,
        verbose=False,
    )

    # Step 6: Compute derived quantities
    derived = compute_derived_quantities(fit, config)

    return TuningForkResult(
        observed=observed,
        fit=fit,
        star_peaks=star_peaks,
        gas_peaks=gas_peaks,
        derived=derived,
        config={
            'lap_min': config.lap_min,
            'lap_max': config.lap_max,
            'n_apertures': config.n_apertures,
            'tstar': config.tstar,
            'peak_prof': config.peak_prof,
            'n_star_peaks': len(star_peaks),
            'n_gas_peaks': len(gas_peaks),
        },
        star_nn=star_nn,
        gas_nn=gas_nn,
    )


def compute_derived_quantities(
    fit: FitResult,
    config: TuningForkConfig,
    tdepl: float = 2000.0,
    fcl: float = 1.0,
    fgmc: float = 1.0,
) -> Dict[str, float]:
    """
    Compute derived physical quantities from fit results.

    Args:
        fit: FitResult from KL14 uncertainty principle fitting
        config: TuningForkConfig used for analysis
        tdepl: Depletion time (Myr) for efficiency calculation
        fcl: Fraction of stellar emission in compact regions
        fgmc: Fraction of gas in giant molecular clouds

    Returns:
        Dictionary of derived quantities
    """
    tgas = fit.tgas
    tover = fit.tover
    lambda_ = fit.lambda_
    tstar = config.tstar

    # Total timeline
    ttotal = tgas + tstar - tover

    # Star formation efficiency
    esf = f_esf(tgas, tdepl, fcl, fgmc)

    # Feedback velocity
    vfb = f_vfb(tover, lambda_)

    # Mass loading factors
    etainst = f_etainst(tgas, tover, esf)
    etaavg = f_etaavg(esf)

    return {
        'tgas': tgas,
        'tover': tover,
        'lambda': lambda_,
        'tstar': tstar,
        'ttotal': ttotal,
        'esf': esf,
        'vfb': vfb,
        'etainst': etainst,
        'etaavg': etaavg,
    }


def generate_model_curve(
    fit: FitResult,
    apertures: Optional[np.ndarray] = None,
    n_points: int = 100,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate smooth model curves for plotting.

    Args:
        fit: FitResult from KL14 uncertainty principle fitting
        apertures: Aperture values to evaluate at (optional)
        n_points: Number of points if apertures not provided

    Returns:
        Tuple of (apertures, model_star, model_gas)
    """
    from heisenberg.core.model import f_fluxratiostar, f_fluxratiogas

    if apertures is None:
        # Generate smooth curve
        lap_min = fit.lambda_ / 5
        lap_max = fit.lambda_ * 10
        apertures = np.logspace(np.log10(lap_min), np.log10(lap_max), n_points)

    tgas = fit.tgas
    tstar = fit.tstar
    tover = fit.tover
    lambda_ = fit.lambda_

    # Use default contrast values
    surfcontrasts = 5.0
    surfcontrastg = 5.0
    beta_star = 1.0
    beta_gas = 1.0
    peak_prof = 2

    model_star = f_fluxratiostar(
        tgas=tgas,
        tstar=tstar,
        tover=tover,
        laps=apertures,
        lambda_=lambda_,
        beta_gas=beta_gas,
        surfcontrasts=surfcontrasts,
        surfcontrastg=surfcontrastg,
        peak_prof=peak_prof,
    )

    model_gas = f_fluxratiogas(
        tgas=tgas,
        tstar=tstar,
        tover=tover,
        lapg=apertures,
        lambda_=lambda_,
        beta_star=beta_star,
        surfcontrasts=surfcontrasts,
        surfcontrastg=surfcontrastg,
        peak_prof=peak_prof,
    )

    return apertures, model_star, model_gas


def generate_plots(
    result: TuningForkResult,
    star_image: np.ndarray,
    gas_image: np.ndarray,
    output_dir: Union[str, Path],
    galaxy_name: str = 'galaxy',
    star_sensitivity: Optional[np.ndarray] = None,
    gas_sensitivity: Optional[np.ndarray] = None,
    dpi: int = 150,
    formats: List[str] = None,
) -> List[Path]:
    """
    Generate all diagnostic plots for a tuning fork analysis.

    This is a convenience wrapper around the plotting module functions.
    Generates tuning fork diagram, residuals, PDFs, maps with peaks,
    and sensitivity histograms.

    Args:
        result: TuningForkResult from run_tuningfork
        star_image: Stellar tracer image
        gas_image: Gas tracer image
        output_dir: Directory for output plots
        galaxy_name: Name for file prefixes
        star_sensitivity: Stellar sensitivity map (optional)
        gas_sensitivity: Gas sensitivity map (optional)
        dpi: Output resolution
        formats: Output file formats (default: ['png', 'pdf'])

    Returns:
        List of generated file paths
    """
    from heisenberg.plotting import generate_all_plots

    if formats is None:
        formats = ['png', 'pdf']

    return generate_all_plots(
        result=result,
        star_image=star_image,
        gas_image=gas_image,
        output_dir=output_dir,
        galaxy_name=galaxy_name,
        star_sensitivity=star_sensitivity,
        gas_sensitivity=gas_sensitivity,
        dpi=dpi,
        formats=formats,
    )


def export_ds9_regions(
    result: TuningForkResult,
    output_dir: Union[str, Path],
    galaxy_name: str = 'galaxy',
) -> List[Path]:
    """
    Export peak positions to DS9 region files.

    Args:
        result: TuningForkResult from run_tuningfork
        output_dir: Directory for output files
        galaxy_name: Name for file prefixes

    Returns:
        List of generated file paths
    """
    from heisenberg.plotting import export_peaks_to_ds9

    return export_peaks_to_ds9(
        result=result,
        output_dir=output_dir,
        galaxy_name=galaxy_name,
    )

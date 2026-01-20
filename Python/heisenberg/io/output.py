"""
Output file generation for Heisenberg analysis results.

Supports writing and reading analysis results in JSON and YAML formats,
preserving all fit parameters, derived quantities, and metadata.

References:
    IDL: write_heisenberg_output.pro
"""

import json
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any, List, Union
import numpy as np

from heisenberg.core.tuningfork import (
    TuningForkResult,
    TuningForkData,
    TuningForkConfig,
)
from heisenberg.core.fitting import FitResult
from heisenberg.peaks import DetectedPeak


@dataclass
class OutputConfig:
    """
    Configuration for output file generation.

    Attributes:
        format: Output format ('json', 'yaml')
        include_arrays: Include probability distribution arrays
        include_metadata: Include analysis metadata (timestamp, version)
        include_peaks: Include individual peak data
        precision: Decimal precision for floating point values
    """
    format: str = "json"
    include_arrays: bool = True
    include_metadata: bool = True
    include_peaks: bool = True
    precision: int = 6


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy arrays and types."""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)


def _fit_result_to_dict(
    fit: FitResult,
    include_arrays: bool = True,
) -> Dict[str, Any]:
    """
    Convert FitResult to dictionary.

    Args:
        fit: FitResult from KL14 uncertainty principle fitting
        include_arrays: Whether to include probability arrays

    Returns:
        Dictionary representation
    """
    result = {
        # Best-fit parameters
        "tgas": float(fit.tgas),
        "tover": float(fit.tover),
        "lambda": float(fit.lambda_),
        "tstar": float(fit.tstar),
        # Asymmetric errors
        "tgas_errmin": float(fit.tgas_errmin),
        "tgas_errmax": float(fit.tgas_errmax),
        "tover_errmin": float(fit.tover_errmin),
        "tover_errmax": float(fit.tover_errmax),
        "lambda_errmin": float(fit.lambda_errmin),
        "lambda_errmax": float(fit.lambda_errmax),
        # Fit quality
        "chi2_min": float(fit.chi2_min),
        "ndof": int(fit.ndof),
        # Auxiliary parameters
        "beta_star": float(fit.beta_star),
        "beta_gas": float(fit.beta_gas),
        "surfcon_star": float(fit.surfcon_star),
        "surfcon_gas": float(fit.surfcon_gas),
        # Correlations
        "corr_tgas_tover": float(fit.corr_tgas_tover),
        "corr_tgas_lambda": float(fit.corr_tgas_lambda),
        "corr_tover_lambda": float(fit.corr_tover_lambda),
    }

    if include_arrays:
        result["tgas_arr"] = fit.tgas_arr.tolist()
        result["tover_arr"] = fit.tover_arr.tolist()
        result["lambda_arr"] = fit.lambda_arr.tolist()
        result["prob_tgas"] = fit.prob_tgas.tolist()
        result["prob_tover"] = fit.prob_tover.tolist()
        result["prob_lambda"] = fit.prob_lambda.tolist()

    return result


def _tuningfork_data_to_dict(data: TuningForkData) -> Dict[str, Any]:
    """Convert TuningForkData to dictionary."""
    return {
        "apertures": data.apertures.tolist(),
        "fluxratio_star": data.fluxratio_star.tolist(),
        "fluxratio_gas": data.fluxratio_gas.tolist(),
        "err_star": data.err_star.tolist(),
        "err_gas": data.err_gas.tolist(),
        "n_star_peaks": data.n_star_peaks.tolist(),
        "n_gas_peaks": data.n_gas_peaks.tolist(),
    }


def _peak_to_dict(peak: DetectedPeak) -> Dict[str, Any]:
    """Convert DetectedPeak to dictionary."""
    result = {
        "x": float(peak.x),
        "y": float(peak.y),
        "total_flux": float(peak.total_flux),
        "npix": int(peak.npix),
    }
    if peak.stats is not None:
        result["stats"] = {
            "peak_flux": float(peak.stats.peak_flux),
            "fwhm_x": float(peak.stats.fwhm_x),
            "fwhm_y": float(peak.stats.fwhm_y),
            "radius": float(peak.stats.radius),
            "on_edge": bool(peak.stats.on_edge),
        }
    return result


def result_to_dict(
    result: TuningForkResult,
    config: Optional[OutputConfig] = None,
) -> Dict[str, Any]:
    """
    Convert TuningForkResult to dictionary for serialization.

    Args:
        result: TuningForkResult from analysis
        config: Output configuration

    Returns:
        Dictionary suitable for JSON/YAML serialization
    """
    if config is None:
        config = OutputConfig()

    output: Dict[str, Any] = {}

    # Metadata
    if config.include_metadata:
        output["metadata"] = {
            "timestamp": datetime.now().isoformat(),
            "version": "1.0.0",
            "format_version": "1",
        }

    # Fit results
    output["fit"] = _fit_result_to_dict(result.fit, config.include_arrays)

    # Observed data
    output["observed"] = _tuningfork_data_to_dict(result.observed)

    # Derived quantities
    output["derived"] = {k: float(v) for k, v in result.derived.items()}

    # Configuration
    output["config"] = result.config

    # Peak data
    if config.include_peaks:
        output["star_peaks"] = [_peak_to_dict(p) for p in result.star_peaks]
        output["gas_peaks"] = [_peak_to_dict(p) for p in result.gas_peaks]

    return output


def write_json(
    result: TuningForkResult,
    path: Union[str, Path],
    config: Optional[OutputConfig] = None,
) -> None:
    """
    Write results to JSON file.

    Args:
        result: TuningForkResult from analysis
        path: Output file path
        config: Output configuration
    """
    if config is None:
        config = OutputConfig(format="json")

    data = result_to_dict(result, config)

    path = Path(path)
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, cls=NumpyEncoder)


def write_yaml(
    result: TuningForkResult,
    path: Union[str, Path],
    config: Optional[OutputConfig] = None,
) -> None:
    """
    Write results to YAML file.

    Args:
        result: TuningForkResult from analysis
        path: Output file path
        config: Output configuration

    Raises:
        ImportError: If PyYAML is not installed
    """
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML required for YAML output: pip install pyyaml")

    if config is None:
        config = OutputConfig(format="yaml")

    data = result_to_dict(result, config)

    path = Path(path)
    with open(path, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def write_results(
    result: TuningForkResult,
    path: Union[str, Path],
    config: Optional[OutputConfig] = None,
) -> None:
    """
    Write analysis results to file.

    Automatically selects format based on file extension or config.

    Args:
        result: TuningForkResult from analysis
        path: Output file path
        config: Output configuration (format auto-detected if not specified)
    """
    path = Path(path)

    if config is None:
        config = OutputConfig()
        # Auto-detect format from extension
        if path.suffix.lower() in ('.yaml', '.yml'):
            config.format = 'yaml'
        else:
            config.format = 'json'

    if config.format == 'yaml':
        write_yaml(result, path, config)
    else:
        write_json(result, path, config)


def _dict_to_fit_result(data: Dict[str, Any]) -> FitResult:
    """Reconstruct FitResult from dictionary."""
    return FitResult(
        tgas=data["tgas"],
        tover=data["tover"],
        lambda_=data["lambda"],
        tstar=data["tstar"],
        tgas_errmin=data["tgas_errmin"],
        tgas_errmax=data["tgas_errmax"],
        tover_errmin=data["tover_errmin"],
        tover_errmax=data["tover_errmax"],
        lambda_errmin=data["lambda_errmin"],
        lambda_errmax=data["lambda_errmax"],
        chi2_min=data["chi2_min"],
        ndof=data["ndof"],
        beta_star=data["beta_star"],
        beta_gas=data["beta_gas"],
        surfcon_star=data["surfcon_star"],
        surfcon_gas=data["surfcon_gas"],
        tgas_arr=np.array(data.get("tgas_arr", [])),
        tover_arr=np.array(data.get("tover_arr", [])),
        lambda_arr=np.array(data.get("lambda_arr", [])),
        prob_tgas=np.array(data.get("prob_tgas", [])),
        prob_tover=np.array(data.get("prob_tover", [])),
        prob_lambda=np.array(data.get("prob_lambda", [])),
        corr_tgas_tover=data["corr_tgas_tover"],
        corr_tgas_lambda=data["corr_tgas_lambda"],
        corr_tover_lambda=data["corr_tover_lambda"],
    )


def _dict_to_tuningfork_data(data: Dict[str, Any]) -> TuningForkData:
    """Reconstruct TuningForkData from dictionary."""
    return TuningForkData(
        apertures=np.array(data["apertures"]),
        fluxratio_star=np.array(data["fluxratio_star"]),
        fluxratio_gas=np.array(data["fluxratio_gas"]),
        err_star=np.array(data["err_star"]),
        err_gas=np.array(data["err_gas"]),
        n_star_peaks=np.array(data["n_star_peaks"]),
        n_gas_peaks=np.array(data["n_gas_peaks"]),
    )


def _dict_to_peak(data: Dict[str, Any]) -> DetectedPeak:
    """Reconstruct DetectedPeak from dictionary."""
    from heisenberg.peaks.statistics import PeakStatistics

    stats = None
    if "stats" in data:
        s = data["stats"]
        stats = PeakStatistics(
            id=0,
            x=data["x"],
            y=data["y"],
            peak_flux=s["peak_flux"],
            total_flux=data["total_flux"],
            fwhm_x=s["fwhm_x"],
            fwhm_y=s["fwhm_y"],
            radius=s["radius"],
            npix=data["npix"],
            on_edge=s["on_edge"],
        )

    return DetectedPeak(
        x=data["x"],
        y=data["y"],
        total_flux=data["total_flux"],
        npix=data["npix"],
        stats=stats,
    )


def read_json(path: Union[str, Path]) -> TuningForkResult:
    """
    Read results from JSON file.

    Args:
        path: Input file path

    Returns:
        TuningForkResult reconstructed from file
    """
    path = Path(path)
    with open(path, 'r') as f:
        data = json.load(f)

    return _dict_to_result(data)


def read_yaml(path: Union[str, Path]) -> TuningForkResult:
    """
    Read results from YAML file.

    Args:
        path: Input file path

    Returns:
        TuningForkResult reconstructed from file

    Raises:
        ImportError: If PyYAML is not installed
    """
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML required for YAML input: pip install pyyaml")

    path = Path(path)
    with open(path, 'r') as f:
        data = yaml.safe_load(f)

    return _dict_to_result(data)


def _dict_to_result(data: Dict[str, Any]) -> TuningForkResult:
    """Reconstruct TuningForkResult from dictionary."""
    fit = _dict_to_fit_result(data["fit"])
    observed = _dict_to_tuningfork_data(data["observed"])

    star_peaks = []
    gas_peaks = []
    if "star_peaks" in data:
        star_peaks = [_dict_to_peak(p) for p in data["star_peaks"]]
    if "gas_peaks" in data:
        gas_peaks = [_dict_to_peak(p) for p in data["gas_peaks"]]

    return TuningForkResult(
        observed=observed,
        fit=fit,
        star_peaks=star_peaks,
        gas_peaks=gas_peaks,
        derived=data.get("derived", {}),
        config=data.get("config", {}),
    )


def read_results(path: Union[str, Path]) -> TuningForkResult:
    """
    Read analysis results from file.

    Automatically selects format based on file extension.

    Args:
        path: Input file path

    Returns:
        TuningForkResult reconstructed from file
    """
    path = Path(path)

    if path.suffix.lower() in ('.yaml', '.yml'):
        return read_yaml(path)
    else:
        return read_json(path)


def write_summary(
    result: TuningForkResult,
    path: Union[str, Path],
) -> None:
    """
    Write a human-readable text summary of results.

    Args:
        result: TuningForkResult from analysis
        path: Output file path
    """
    path = Path(path)
    fit = result.fit

    lines = [
        "=" * 60,
        "HEISENBERG ANALYSIS RESULTS",
        "=" * 60,
        "",
        "Best-fit parameters:",
        f"  t_gas    = {fit.tgas:.2f} (+{fit.tgas_errmax:.2f}/-{fit.tgas_errmin:.2f}) Myr",
        f"  t_over   = {fit.tover:.2f} (+{fit.tover_errmax:.2f}/-{fit.tover_errmin:.2f}) Myr",
        f"  lambda   = {fit.lambda_:.1f} (+{fit.lambda_errmax:.1f}/-{fit.lambda_errmin:.1f}) pc",
        f"  t_star   = {fit.tstar:.2f} Myr (input)",
        "",
        f"Fit quality: chi2_min = {fit.chi2_min:.2f}, ndof = {fit.ndof}",
        "",
        "Derived quantities:",
    ]

    for key, value in result.derived.items():
        lines.append(f"  {key} = {value:.4g}")

    lines.extend([
        "",
        f"Number of stellar peaks: {result.config.get('n_star_peaks', 'N/A')}",
        f"Number of gas peaks: {result.config.get('n_gas_peaks', 'N/A')}",
        "",
        "=" * 60,
    ])

    with open(path, 'w') as f:
        f.write('\n'.join(lines))

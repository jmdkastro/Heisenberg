"""
Interactive peak finding TUI.

Provides a text user interface for iteratively adjusting peak finding
parameters until the user is satisfied with the identified peaks.

This module translates the IDL interactive_peak_find.pro functionality.

References:
    IDL: interactive_peak_find.pro
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any, Callable, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from heisenberg.peaks.detection import DetectedPeak, PeakDetectionConfig


@dataclass
class InteractivePeakConfig:
    """
    Configuration for interactive peak finding.

    Attributes:
        npixmin: Minimum number of pixels per peak
        nsigma: Detection threshold (multiples of sensitivity)
        loglevels: Use logarithmic contour levels
        logrange_s: Log range for stellar contour levels
        logspacing_s: Log spacing for stellar contour levels
        logrange_g: Log range for gas contour levels
        logspacing_g: Log spacing for gas contour levels
        nlinlevel_s: Number of linear levels for stellar peaks
        nlinlevel_g: Number of linear levels for gas peaks
        flux_weighted: Use flux-weighted positions
    """
    npixmin: int = 20
    nsigma: float = 5.0
    loglevels: bool = True
    logrange_s: float = 2.0
    logspacing_s: float = 0.5
    logrange_g: float = 2.0
    logspacing_g: float = 0.5
    nlinlevel_s: int = 11
    nlinlevel_g: int = 11
    flux_weighted: bool = False


@dataclass
class InteractiveResult:
    """
    Result from interactive peak finding.

    Attributes:
        star_peaks: Detected stellar peaks
        gas_peaks: Detected gas peaks
        config: Final configuration used
        finalized: Whether parameters are finalized for future iterations
        n_iterations: Number of iterations performed
    """
    star_peaks: List['DetectedPeak']
    gas_peaks: List['DetectedPeak']
    config: InteractivePeakConfig
    finalized: bool = False
    n_iterations: int = 0


def _print_current_settings(config: InteractivePeakConfig, iteration: int) -> None:
    """Print current peak finding settings."""
    print(f"\nPeakID #{iteration}. Current settings:")
    print(f"  npixmin          {config.npixmin}")
    print(f"  nsigma           {config.nsigma}")
    print(f"  logrange_s       {config.logrange_s}")
    print(f"  logspacing_s     {config.logspacing_s}")
    print(f"  logrange_g       {config.logrange_g}")
    print(f"  logspacing_g     {config.logspacing_g}")
    print(f"  nlinlevel_s      {config.nlinlevel_s}")
    print(f"  nlinlevel_g      {config.nlinlevel_g}")


def _prompt_yes_no(prompt: str) -> bool:
    """Prompt user for yes/no response."""
    while True:
        response = input(f"{prompt} (yes/no): ").strip().lower()
        if response in ('y', 'yes'):
            return True
        elif response in ('n', 'no'):
            return False
        else:
            print("Please enter 'yes' or 'no'")


def _prompt_value(name: str, current_value: float, is_global: bool = False) -> float:
    """
    Prompt user to change a parameter value.

    Args:
        name: Parameter name
        current_value: Current value
        is_global: Whether this is a global parameter affecting both star/gas

    Returns:
        New value (or current if unchanged)
    """
    global_note = " (affects both star and gas)" if is_global else ""
    prompt = f"  {name}{global_note} [{current_value}]: "
    response = input(prompt).strip()

    if not response:
        return current_value

    try:
        return float(response)
    except ValueError:
        print(f"  Invalid value, keeping {current_value}")
        return current_value


def _prompt_int_value(name: str, current_value: int, is_global: bool = False) -> int:
    """Prompt user to change an integer parameter value."""
    global_note = " (affects both star and gas)" if is_global else ""
    prompt = f"  {name}{global_note} [{current_value}]: "
    response = input(prompt).strip()

    if not response:
        return current_value

    try:
        return int(response)
    except ValueError:
        print(f"  Invalid value, keeping {current_value}")
        return current_value


def _run_peak_detection(
    image: np.ndarray,
    config: InteractivePeakConfig,
    sensitivity: Optional[np.ndarray] = None,
    is_star: bool = True,
) -> List['DetectedPeak']:
    """
    Run peak detection with current configuration.

    Args:
        image: Image to find peaks in
        config: Current configuration
        sensitivity: Sensitivity map
        is_star: True for stellar map, False for gas map

    Returns:
        List of detected peaks
    """
    from heisenberg.peaks.detection import PeakDetectionConfig, find_peaks

    if config.loglevels:
        logrange = config.logrange_s if is_star else config.logrange_g
        logspacing = config.logspacing_s if is_star else config.logspacing_g
        # Calculate nlevels like IDL tuningfork.pro: nlevels = logrange / logspacing + 1
        nlevels = int(logrange / logspacing) + 1
    else:
        logrange = 2.0
        logspacing = 0.5
        nlevels = config.nlinlevel_s if is_star else config.nlinlevel_g

    peak_config = PeakDetectionConfig(
        npixmin=config.npixmin,
        nsigma=config.nsigma,
        loglevels=config.loglevels,
        logrange=logrange,
        logspacing=logspacing,
        nlevels=nlevels,
        flux_weighted=config.flux_weighted,
    )

    return find_peaks(image, config=peak_config, sensitivity=sensitivity)


def _display_peaks_matplotlib(
    image: np.ndarray,
    peaks: List['DetectedPeak'],
    title: str,
    ax=None,
) -> None:
    """Display image with peaks using matplotlib."""
    try:
        import matplotlib.pyplot as plt
        from heisenberg.plotting.image_display import display_image

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))

        display_image(image, ax=ax, stretch='log')

        if peaks:
            x_coords = [p.x for p in peaks]
            y_coords = [p.y for p in peaks]
            ax.scatter(x_coords, y_coords, marker='x', s=50, c='red', linewidths=1.5)

        ax.set_title(f"{title} ({len(peaks)} peaks)")
        plt.show(block=False)
        plt.pause(0.1)

    except ImportError:
        print(f"[matplotlib not available] {title}: {len(peaks)} peaks")


def _display_peaks_ds9(
    image_path: Path,
    peaks: List['DetectedPeak'],
    title: str,
) -> None:
    """Display image with peaks in DS9."""
    try:
        import subprocess

        # Check if DS9 is available
        result = subprocess.run(['which', 'ds9'], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[DS9 not found] {title}: {len(peaks)} peaks")
            return

        # Launch DS9 with the image
        subprocess.Popen(
            ['ds9', str(image_path), '-title', title, '-scale', 'log'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Export peaks to region file if there are any
        if peaks:
            from heisenberg.plotting.ds9_regions import write_peak_regions
            region_path = Path(f"/tmp/{title.replace(' ', '_')}.reg")
            write_peak_regions(peaks, region_path, color='red')
            print(f"Region file: {region_path}")

    except Exception as e:
        print(f"[DS9 error: {e}] {title}: {len(peaks)} peaks")


def _export_region_file(
    peaks: List['DetectedPeak'],
    output_path: Path,
    color: str = 'red',
) -> None:
    """Export peaks to DS9 region file."""
    from heisenberg.plotting.ds9_regions import write_peak_regions
    write_peak_regions(peaks, output_path, color=color)


def _write_report(
    config: InteractivePeakConfig,
    output_dir: Path,
    galaxy_name: str,
    finalized: bool,
) -> Path:
    """Write interactive peak finding report."""
    report_path = output_dir / f"{galaxy_name}_interactive_peak_ID_report.dat"

    with open(report_path, 'w') as f:
        f.write("# Interactive Peak Identification report. Final selected parameters\n")
        f.write(f"npixmin          {config.npixmin}\n")
        f.write(f"nsigma           {config.nsigma}\n")
        f.write(f"logrange_s       {config.logrange_s}\n")
        f.write(f"logspacing_s     {config.logspacing_s}\n")
        f.write(f"logrange_g       {config.logrange_g}\n")
        f.write(f"logspacing_g     {config.logspacing_g}\n")
        f.write(f"nlinlevel_s      {config.nlinlevel_s}\n")
        f.write(f"nlinlevel_g      {config.nlinlevel_g}\n")
        f.write(f"finalised        {1 if finalized else 0}\n")

    return report_path


def interactive_peak_find(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    config: Optional[InteractivePeakConfig] = None,
    star_sensitivity: Optional[np.ndarray] = None,
    gas_sensitivity: Optional[np.ndarray] = None,
    output_dir: Optional[Path] = None,
    galaxy_name: str = 'galaxy',
    use_ds9: bool = False,
    use_matplotlib: bool = True,
) -> InteractiveResult:
    """
    Run interactive peak finding TUI.

    This provides a text user interface for iteratively adjusting peak
    finding parameters until the user is satisfied with the results.

    Args:
        star_image: Stellar tracer image
        gas_image: Gas tracer image
        config: Initial configuration (uses defaults if None)
        star_sensitivity: Stellar sensitivity map
        gas_sensitivity: Gas sensitivity map
        output_dir: Output directory for region files and report
        galaxy_name: Galaxy name for output files
        use_ds9: Display images in DS9
        use_matplotlib: Display images with matplotlib

    Returns:
        InteractiveResult with final peaks and configuration
    """
    if config is None:
        config = InteractivePeakConfig()

    if output_dir is None:
        output_dir = Path(".")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    iteration = 0
    run_stars = True
    run_gas = True
    star_peaks: List['DetectedPeak'] = []
    gas_peaks: List['DetectedPeak'] = []

    print("\n" + "=" * 60)
    print("Interactive Peak Finding")
    print("=" * 60)

    while True:
        # Run peak detection
        if run_stars:
            star_peaks = _run_peak_detection(
                star_image, config, star_sensitivity, is_star=True
            )

        if run_gas:
            gas_peaks = _run_peak_detection(
                gas_image, config, gas_sensitivity, is_star=False
            )

        # Check for too few peaks
        empty_stars = len(star_peaks) < 2
        empty_gas = len(gas_peaks) < 2

        # Print current settings
        _print_current_settings(config, iteration)

        # Display peaks
        print(f"\nStar peaks: {len(star_peaks)}")
        print(f"Gas peaks: {len(gas_peaks)}")

        if use_matplotlib:
            try:
                import matplotlib.pyplot as plt
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
                _display_peaks_matplotlib(star_image, star_peaks, f"Stellar #{iteration}", ax1)
                _display_peaks_matplotlib(gas_image, gas_peaks, f"Gas #{iteration}", ax2)
                plt.tight_layout()
                plt.show(block=False)
                plt.pause(0.1)
            except ImportError:
                pass

        # Check if too few peaks
        if empty_stars or empty_gas:
            if empty_stars and empty_gas:
                empty_msg = "star and gas maps"
            elif empty_stars:
                empty_msg = "star map"
            else:
                empty_msg = "gas map"
            print(f"\nToo few peaks found in {empty_msg}. Must adjust parameters.")
            satisfied = False
        else:
            # Ask if satisfied
            print("")
            satisfied = _prompt_yes_no("Are you satisfied with the identified peaks?")

        if satisfied:
            break

        # Adjust parameters
        print("\nAdjust parameters (press Enter to keep current value):")

        old_config = InteractivePeakConfig(
            npixmin=config.npixmin,
            nsigma=config.nsigma,
            logrange_s=config.logrange_s,
            logspacing_s=config.logspacing_s,
            logrange_g=config.logrange_g,
            logspacing_g=config.logspacing_g,
            nlinlevel_s=config.nlinlevel_s,
            nlinlevel_g=config.nlinlevel_g,
        )

        # Global parameters
        config.npixmin = _prompt_int_value("npixmin", config.npixmin, is_global=True)
        config.nsigma = _prompt_value("nsigma", config.nsigma, is_global=True)

        # Determine what to rerun
        run_stars = False
        run_gas = False

        if config.npixmin != old_config.npixmin or config.nsigma != old_config.nsigma:
            run_stars = True
            run_gas = True

        if config.loglevels:
            # Log spacing parameters
            config.logrange_s = _prompt_value("logrange_s", config.logrange_s)
            config.logspacing_s = _prompt_value("logspacing_s", config.logspacing_s)
            if config.logrange_s != old_config.logrange_s or config.logspacing_s != old_config.logspacing_s:
                run_stars = True

            config.logrange_g = _prompt_value("logrange_g", config.logrange_g)
            config.logspacing_g = _prompt_value("logspacing_g", config.logspacing_g)
            if config.logrange_g != old_config.logrange_g or config.logspacing_g != old_config.logspacing_g:
                run_gas = True
        else:
            # Linear level parameters
            config.nlinlevel_s = _prompt_int_value("nlinlevel_s", config.nlinlevel_s)
            if config.nlinlevel_s != old_config.nlinlevel_s:
                run_stars = True

            config.nlinlevel_g = _prompt_int_value("nlinlevel_g", config.nlinlevel_g)
            if config.nlinlevel_g != old_config.nlinlevel_g:
                run_gas = True

        if not run_stars and not run_gas:
            print("\nNo parameters changed. Please change at least one parameter to rerun.")
            continue

        iteration += 1

    # Ask about finalizing
    print("")
    finalized = _prompt_yes_no(
        "Do you want to finalize peak finding parameters for future iterations?\n"
        "  'yes' = do not adjust in future, 'no' = continue adjusting"
    )

    if finalized:
        print("Interactive peak finding will be disabled in future iterations.")
    else:
        print("Interactive peak finding will continue in future iterations.")

    # Export region files
    star_region_path = output_dir / f"{galaxy_name}_star_peaks.reg"
    gas_region_path = output_dir / f"{galaxy_name}_gas_peaks.reg"
    _export_region_file(star_peaks, star_region_path, color='red')
    _export_region_file(gas_peaks, gas_region_path, color='cyan')
    print(f"\nRegion files saved to {output_dir}")

    # Write report
    report_path = _write_report(config, output_dir, galaxy_name, finalized)
    print(f"Report saved to {report_path}")

    return InteractiveResult(
        star_peaks=star_peaks,
        gas_peaks=gas_peaks,
        config=config,
        finalized=finalized,
        n_iterations=iteration + 1,
    )


def non_interactive_peak_find(
    star_image: np.ndarray,
    gas_image: np.ndarray,
    config: Optional[InteractivePeakConfig] = None,
    star_sensitivity: Optional[np.ndarray] = None,
    gas_sensitivity: Optional[np.ndarray] = None,
) -> Tuple[List['DetectedPeak'], List['DetectedPeak']]:
    """
    Run non-interactive peak finding.

    This is equivalent to running interactive_peak_find with peak_find_tui=0
    in the IDL version.

    Args:
        star_image: Stellar tracer image
        gas_image: Gas tracer image
        config: Configuration
        star_sensitivity: Stellar sensitivity map
        gas_sensitivity: Gas sensitivity map

    Returns:
        Tuple of (star_peaks, gas_peaks)
    """
    if config is None:
        config = InteractivePeakConfig()

    star_peaks = _run_peak_detection(
        star_image, config, star_sensitivity, is_star=True
    )
    gas_peaks = _run_peak_detection(
        gas_image, config, gas_sensitivity, is_star=False
    )

    return star_peaks, gas_peaks

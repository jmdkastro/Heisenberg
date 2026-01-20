"""
Command-line interface for Heisenberg.

This module provides the CLI for running Heisenberg analyses from the command line.

Usage:
    heisenberg run input_file
    heisenberg run input_file --no-diffuse-filtering
    heisenberg validate input_file
    heisenberg create-input output_file
"""

import click
from pathlib import Path
from typing import Optional


@click.group()
@click.version_option()
def main():
    """
    Heisenberg: Uncertainty Principle for Star Formation.

    A tool for measuring GMC lifetime, gas clearance timescale, and star
    formation efficiency from galaxy tracer maps.
    """
    pass


def _config_to_tuningfork_config(config):
    """Convert HeisenbergConfig to TuningForkConfig."""
    from heisenberg.core.tuningfork import TuningForkConfig

    return TuningForkConfig(
        lap_min=config.aperture.lapmin,
        lap_max=config.aperture.lapmax,
        n_apertures=config.aperture.naperture,
        n_mc=config.fitting.nmc,
        tstar=config.timeline.tstariso,
        peak_prof=config.flags4.peak_prof,
        npixmin_star=config.peak_id.npixmin,
        npixmin_gas=config.peak_id.npixmin,
        nsigma_star=config.peak_id.nsigma,
        nsigma_gas=config.peak_id.nsigma,
        tstar_incl=config.flags4.tstar_incl,
        tgas_min=config.timeline.tgasmini,
        tgas_max=config.timeline.tgasmaxi,
        tover_min=config.timeline.tovermini,
        ndepth=config.fitting.ndepth,
        ntry=config.fitting.ntry,
    )


def _config_to_diffuse_config(config):
    """Convert HeisenbergConfig to DiffuseIterationConfig."""
    from heisenberg.core.diffuse_iteration import DiffuseIterationConfig

    return DiffuseIterationConfig(
        use_guess=config.iteration.use_guess,
        initial_guess=config.iteration.initial_guess,
        iter_criterion=config.iteration.iter_criterion,
        iter_crit_len=config.iteration.iter_crit_len,
        iter_nmax=config.iteration.iter_nmax,
        iter_filter=config.iteration.iter_filter,
        iter_bwo=config.iteration.iter_bwo,
        iter_len_conv=config.iteration.iter_len_conv,
        use_noisecut=config.noise.use_noisecut,
        noisethresh_s=config.noise.noisethresh_s,
        noisethresh_g=config.noise.noisethresh_g,
    )


def _load_data_from_config(config, verbose: bool = False):
    """Load FITS data based on configuration."""
    from heisenberg.io.fits import read_fits, get_pixel_scale_pc

    datadir = config.files.datadir

    # Load stellar map
    star_path = datadir / config.files.starfile
    if verbose:
        click.echo(f"Loading stellar map: {star_path}")
    star_data, star_wcs, _ = read_fits(star_path)

    # Load gas map
    gas_path = datadir / config.files.gasfile
    if verbose:
        click.echo(f"Loading gas map: {gas_path}")
    gas_data, gas_wcs, _ = read_fits(gas_path)

    # Get pixel scale in pc
    if star_wcs is not None:
        pixel_scale_pc = get_pixel_scale_pc(star_wcs, config.basic_map.distance)
    else:
        raise ValueError("No WCS found - cannot determine pixel scale")

    if verbose:
        click.echo(f"Image size: {star_data.shape}")
        click.echo(f"Pixel scale: {pixel_scale_pc:.2f} pc/pixel")

    return star_data, gas_data, pixel_scale_pc


def _print_results(result, verbose: bool = False):
    """Print analysis results."""
    click.echo("")
    click.echo("=" * 60)
    click.echo("RESULTS")
    click.echo("=" * 60)
    click.echo(f"t_gas    = {result.fit.tgas:.2f} (+{result.fit.tgas_errmax:.2f}/-{result.fit.tgas_errmin:.2f}) Myr")
    click.echo(f"t_over   = {result.fit.tover:.2f} (+{result.fit.tover_errmax:.2f}/-{result.fit.tover_errmin:.2f}) Myr")
    click.echo(f"lambda   = {result.fit.lambda_:.1f} (+{result.fit.lambda_errmax:.1f}/-{result.fit.lambda_errmin:.1f}) pc")
    click.echo(f"chi2_min = {result.fit.chi2_min:.2f}")
    click.echo(f"Peaks: {len(result.star_peaks)} stellar, {len(result.gas_peaks)} gas")

    if verbose and result.derived:
        click.echo("")
        click.echo("Derived quantities:")
        for key, value in result.derived.items():
            click.echo(f"  {key}: {value:.4g}")


@main.command()
@click.argument('input_file', type=click.Path(exists=True, path_type=Path))
@click.option(
    '--no-diffuse-filtering', '-n',
    is_flag=True,
    help='Skip iterative diffuse filtering (run tuningfork directly)'
)
@click.option(
    '--output', '-o',
    type=click.Path(path_type=Path),
    help='Output file path (default: <galaxy>_results.json)'
)
@click.option(
    '--output-dir', '-d',
    type=click.Path(path_type=Path),
    help='Output directory for plots and region files'
)
@click.option(
    '--plot', '-p',
    is_flag=True,
    help='Generate diagnostic plots'
)
@click.option(
    '--regions', '-r',
    is_flag=True,
    help='Export peak positions to DS9 region files'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Enable verbose output'
)
def run(
    input_file: Path,
    no_diffuse_filtering: bool,
    output: Optional[Path],
    output_dir: Optional[Path],
    plot: bool,
    regions: bool,
    verbose: bool,
):
    """
    Run the Heisenberg analysis.

    INPUT_FILE is the path to the configuration file in IDL input_file format.
    """
    from heisenberg.config import load_config
    from heisenberg.io import write_results, write_summary

    click.echo(f"Loading configuration from {input_file}...")
    config = load_config(input_file)
    click.echo(f"Galaxy: {config.galaxy}")
    click.echo(f"Distance: {config.distance:.0f} pc")

    # Set default output paths
    if output is None:
        output = Path(f"{config.galaxy}_results.json")
    if output_dir is None:
        output_dir = Path(".")

    # Load FITS data
    try:
        star_data, gas_data, pixel_scale_pc = _load_data_from_config(config, verbose)
    except Exception as e:
        click.echo(click.style(f"Error loading data: {e}", fg='red'))
        raise SystemExit(1)

    if no_diffuse_filtering:
        # Single-pass analysis (heisenberg_nodf equivalent)
        click.echo("Running single-pass analysis (no diffuse filtering)...")
        from heisenberg.core.tuningfork import run_tuningfork

        tf_config = _config_to_tuningfork_config(config)
        if verbose:
            click.echo(f"Apertures: {tf_config.n_apertures} ({tf_config.lap_min:.0f} - {tf_config.lap_max:.0f} pc)")
            click.echo(f"t_star: {tf_config.tstar:.1f} Myr")

        try:
            result = run_tuningfork(
                star_image=star_data,
                gas_image=gas_data,
                pixel_scale=pixel_scale_pc,
                config=tf_config,
            )
        except Exception as e:
            click.echo(click.style(f"Analysis error: {e}", fg='red'))
            raise SystemExit(1)

    else:
        # Iterative diffuse filtering analysis (heisenberg equivalent)
        click.echo("Running iterative diffuse filtering analysis...")
        from heisenberg.core.diffuse_iteration import run_diffuse_iteration

        tf_config = _config_to_tuningfork_config(config)
        df_config = _config_to_diffuse_config(config)

        if verbose:
            click.echo(f"Max iterations: {df_config.iter_nmax}")
            click.echo(f"Convergence criterion: {df_config.iter_criterion:.1%}")

        try:
            result, history = run_diffuse_iteration(
                star_image=star_data,
                gas_image=gas_data,
                star_pix_to_pc=pixel_scale_pc,
                gas_pix_to_pc=pixel_scale_pc,
                config=df_config,
                tuningfork_config=tf_config,
                verbose=verbose,
            )
            if verbose and history is not None:
                click.echo(f"Converged after {len(history.lambda_values)} iterations")
        except Exception as e:
            click.echo(click.style(f"Analysis error: {e}", fg='red'))
            raise SystemExit(1)

    # Print results
    _print_results(result, verbose)

    # Save results
    write_results(result, output)
    click.echo(f"Results saved to {output}")

    summary_path = output.with_suffix('.txt')
    write_summary(result, summary_path)
    click.echo(f"Summary saved to {summary_path}")

    # Generate plots if requested
    if plot:
        try:
            from heisenberg.core.tuningfork import generate_plots
            output_dir.mkdir(parents=True, exist_ok=True)
            plot_files = generate_plots(
                result=result,
                star_image=star_data,
                gas_image=gas_data,
                output_dir=output_dir,
                galaxy_name=config.galaxy,
            )
            click.echo(f"Plots saved to {output_dir} ({len(plot_files)} files)")
        except ImportError:
            click.echo("Note: matplotlib not available for plotting")
        except Exception as e:
            click.echo(click.style(f"Plotting error: {e}", fg='yellow'))

    # Export DS9 regions if requested
    if regions:
        try:
            from heisenberg.core.tuningfork import export_ds9_regions
            output_dir.mkdir(parents=True, exist_ok=True)
            region_files = export_ds9_regions(
                result=result,
                output_dir=output_dir,
                galaxy_name=config.galaxy,
            )
            click.echo(f"Region files saved to {output_dir} ({len(region_files)} files)")
        except Exception as e:
            click.echo(click.style(f"Region export error: {e}", fg='yellow'))


@main.command()
@click.argument('input_file', type=click.Path(exists=True, path_type=Path))
def validate(input_file: Path):
    """
    Validate an input configuration file.

    Checks that all required parameters are present and valid.
    """
    from heisenberg.config import load_config

    try:
        config = load_config(input_file)
        click.echo(click.style("Configuration is valid!", fg='green'))
        click.echo(f"  Galaxy: {config.galaxy}")
        click.echo(f"  Distance: {config.distance:.0f} pc")
        click.echo(f"  Data directory: {config.datadir}")
        click.echo(f"  Star file: {config.files.starfile}")
        click.echo(f"  Gas file: {config.files.gasfile}")
        click.echo(f"  Apertures: {config.aperture.naperture} ({config.aperture.lapmin:.0f} - {config.aperture.lapmax:.0f} pc)")
    except Exception as e:
        click.echo(click.style(f"Configuration error: {e}", fg='red'))
        raise SystemExit(1)


@main.command('create-input')
@click.argument('output_file', type=click.Path(path_type=Path))
@click.option(
    '--template', '-t',
    type=click.Choice(['minimal', 'full']),
    default='full',
    help='Template type: minimal (required only) or full (all parameters)'
)
def create_input(output_file: Path, template: str):
    """
    Create a template input configuration file.

    OUTPUT_FILE is the path where the template will be written.
    """
    # TODO: Implement template generation
    click.echo(f"Creating {template} template at {output_file}...")
    click.echo("Template generation not yet implemented.")


@main.command('quick-test')
@click.option(
    '--n-peaks', '-n',
    default=15,
    help='Number of peaks to generate'
)
@click.option(
    '--output', '-o',
    type=click.Path(path_type=Path),
    help='Output file path'
)
@click.option(
    '--plot', '-p',
    is_flag=True,
    help='Generate diagnostic plot'
)
def quick_test(n_peaks: int, output: Path, plot: bool):
    """
    Run a quick test with synthetic data.

    This generates synthetic star and gas maps and runs the full
    tuningfork analysis pipeline for testing.
    """
    import numpy as np
    from heisenberg.core.tuningfork import TuningForkConfig, run_tuningfork
    from heisenberg.io import write_results, write_summary

    click.echo(f"Generating synthetic galaxy with {n_peaks} peaks...")

    # Generate synthetic data
    np.random.seed(42)
    shape = (200, 200)
    margin = 30

    def make_gaussian(shape, cy, cx, amp, sigma):
        y, x = np.ogrid[:shape[0], :shape[1]]
        return amp * np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * sigma**2))

    star_image = np.zeros(shape)
    gas_image = np.zeros(shape)

    for i in range(n_peaks):
        y = np.random.randint(margin, shape[0] - margin)
        x = np.random.randint(margin, shape[1] - margin)
        amp = np.random.uniform(50, 150)
        star_image += make_gaussian(shape, y, x, amp, 8.0)
        gas_image += make_gaussian(
            shape,
            y + np.random.randint(-5, 5),
            x + np.random.randint(-5, 5),
            amp * np.random.uniform(0.8, 1.2),
            8.0
        )

    star_image += np.random.normal(0, 5, shape)
    gas_image += np.random.normal(0, 5, shape)

    click.echo("Running tuningfork analysis...")

    config = TuningForkConfig(
        lap_min=20.0,
        lap_max=100.0,
        n_apertures=6,
        n_mc=20,
        npixmin_star=15,
        npixmin_gas=15,
        nsigma_star=3.0,
        nsigma_gas=3.0,
        seed=42,
    )

    result = run_tuningfork(
        star_image,
        gas_image,
        pixel_scale=1.0,
        config=config,
    )

    click.echo("\n" + "=" * 50)
    click.echo("RESULTS")
    click.echo("=" * 50)
    click.echo(f"t_gas    = {result.fit.tgas:.2f} (+{result.fit.tgas_errmax:.2f}/-{result.fit.tgas_errmin:.2f}) Myr")
    click.echo(f"t_over   = {result.fit.tover:.2f} (+{result.fit.tover_errmax:.2f}/-{result.fit.tover_errmin:.2f}) Myr")
    click.echo(f"lambda   = {result.fit.lambda_:.1f} (+{result.fit.lambda_errmax:.1f}/-{result.fit.lambda_errmin:.1f}) pc")
    click.echo(f"chi2_min = {result.fit.chi2_min:.2f}")
    click.echo(f"Peaks: {len(result.star_peaks)} stellar, {len(result.gas_peaks)} gas")

    if output:
        write_results(result, output)
        click.echo(f"\nResults saved to {output}")
        summary_path = output.with_suffix('.txt')
        write_summary(result, summary_path)
        click.echo(f"Summary saved to {summary_path}")

    if plot:
        try:
            from heisenberg.plotting import plot_summary, save_figure
            fig = plot_summary(result, suptitle='Quick Test Results')
            plot_path = (output or Path('quick_test')).with_suffix('.png')
            save_figure(fig, plot_path)
            click.echo(f"Plot saved to {plot_path}")
        except ImportError:
            click.echo("Note: matplotlib not available for plotting")


if __name__ == '__main__':
    main()

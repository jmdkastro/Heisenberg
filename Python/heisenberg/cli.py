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


@click.group()
@click.version_option()
def main():
    """
    Heisenberg: Uncertainty Principle for Star Formation.

    A tool for measuring GMC lifetime, gas clearance timescale, and star
    formation efficiency from galaxy tracer maps.
    """
    pass


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
    help='Output file path (default: galaxy_results.json)'
)
@click.option(
    '--plot', '-p',
    is_flag=True,
    help='Generate diagnostic plots'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Enable verbose output'
)
def run(
    input_file: Path,
    no_diffuse_filtering: bool,
    output: Path,
    plot: bool,
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

    if no_diffuse_filtering:
        click.echo("Running single-pass analysis (no diffuse filtering)...")
        # TODO: Implement when FITS I/O is complete
        click.echo("Note: Full pipeline requires FITS files - not yet implemented.")
        click.echo("Use --help for available options.")
    else:
        click.echo("Running with iterative diffuse filtering...")
        # TODO: Implement iterative analysis (diffuse_iteration equivalent)
        click.echo("Note: Iterative analysis not yet implemented.")
        click.echo("Use -n flag for single-pass analysis.")


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

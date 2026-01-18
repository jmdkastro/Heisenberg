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
    '--verbose', '-v',
    is_flag=True,
    help='Enable verbose output'
)
def run(input_file: Path, no_diffuse_filtering: bool, verbose: bool):
    """
    Run the Heisenberg analysis.

    INPUT_FILE is the path to the configuration file in IDL input_file format.
    """
    from heisenberg.config import load_config

    click.echo(f"Loading configuration from {input_file}...")
    config = load_config(input_file)
    click.echo(f"Galaxy: {config.galaxy}")
    click.echo(f"Distance: {config.distance:.0f} pc")

    if no_diffuse_filtering:
        click.echo("Running without iterative diffuse filtering...")
        # TODO: Implement single-pass analysis (tuningfork equivalent)
        click.echo("Pipeline not yet implemented - foundation only.")
    else:
        click.echo("Running with iterative diffuse filtering...")
        # TODO: Implement iterative analysis (diffuse_iteration equivalent)
        click.echo("Pipeline not yet implemented - foundation only.")


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


if __name__ == '__main__':
    main()

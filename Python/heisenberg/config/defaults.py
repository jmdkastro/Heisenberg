"""
Default parameter values for Heisenberg configuration.

This module provides factory functions to create configuration objects
with sensible default values, making it easy to create configs programmatically.
"""

from pathlib import Path
from typing import Optional

from heisenberg.config.parameters import (
    HeisenbergConfig,
    FileNames,
    PeakIdFileNames,
    MaskFileNames,
    UnfilteredFileNames,
    Flags1,
    Flags2,
    Flags3,
    Flags4,
    BasicMapParams,
    ApertureParams,
    PeakIdParams,
    TimelineParams,
    FittingParams,
    FourierFilterParams,
    ConversionParams,
    SensitivityParams,
    NoiseThresholdParams,
    IterationParams,
)


def create_minimal_config(
    datadir: str | Path,
    galaxy: str,
    starfile: str,
    gasfile: str,
    distance: float,
    centrex: int,
    centrey: int,
    maxradius: float,
) -> HeisenbergConfig:
    """
    Create a minimal configuration with required parameters only.

    This is useful for simple analyses where most defaults are acceptable.

    Args:
        datadir: Path to the directory containing FITS files
        galaxy: Name for this analysis (used in output filenames)
        starfile: Name of the stellar/SF map FITS file
        gasfile: Name of the gas map FITS file
        distance: Distance to the galaxy in parsecs
        centrex: X pixel coordinate of galaxy center (0-indexed)
        centrey: Y pixel coordinate of galaxy center (0-indexed)
        maxradius: Maximum radius for analysis in parsecs

    Returns:
        HeisenbergConfig with minimal required parameters

    Example:
        >>> config = create_minimal_config(
        ...     datadir="/data/ngc300",
        ...     galaxy="NGC300",
        ...     starfile="halpha.fits",
        ...     gasfile="co.fits",
        ...     distance=2.0e6,  # 2 Mpc
        ...     centrex=512,
        ...     centrey=512,
        ...     maxradius=5000.0,  # 5 kpc
        ... )
    """
    return HeisenbergConfig(
        files=FileNames(
            datadir=Path(datadir),
            galaxy=galaxy,
            starfile=starfile,
            gasfile=gasfile,
        ),
        basic_map=BasicMapParams(
            distance=distance,
            centrex=centrex,
            centrey=centrey,
            maxradius=maxradius,
        ),
    )


def create_ngc300_template() -> HeisenbergConfig:
    """
    Create a configuration template based on NGC300 analysis.

    This provides a starting point for nearby galaxy analyses
    with parameters similar to those used in Kruijssen+2019.

    Returns:
        HeisenbergConfig template for NGC300-like galaxies
    """
    return HeisenbergConfig(
        files=FileNames(
            datadir=Path("/path/to/ngc300/data"),
            galaxy="NGC300",
            starfile="halpha.fits",
            gasfile="co21.fits",
        ),
        basic_map=BasicMapParams(
            distance=2.0e6,  # 2 Mpc in pc
            inclination=42.0,
            posangle=111.0,
            centrex=512,
            centrey=512,
            minradius=0.0,
            maxradius=5000.0,  # 5 kpc
        ),
        aperture=ApertureParams(
            lapmin=50.0,
            lapmax=1600.0,
            naperture=7,
            peak_res=0,
            max_res=6,
        ),
        timeline=TimelineParams(
            tstariso=4.0,  # H-alpha reference timescale
            tstariso_errmin=0.5,
            tstariso_errmax=0.5,
        ),
        conversion=ConversionParams(
            convstar=-3.7,  # H-alpha to SFR
            convgas=2.3,    # CO to gas mass
        ),
    )


# YAML template for easy editing
YAML_TEMPLATE = '''# Heisenberg Configuration File
# ================================
# This YAML file configures the Heisenberg analysis.
# Edit the values below to match your data.

# File locations
files:
  datadir: /path/to/data        # Directory containing FITS files
  galaxy: MyGalaxy              # Name for output files
  starfile: star_map.fits       # Stellar/SF tracer map
  gasfile: gas_map.fits         # Gas tracer map
  # Optional secondary maps (uncomment if needed):
  # starfile2: star_peakid.fits
  # gasfile2: gas_peakid.fits
  # starfile3: star_ratio.fits

# Galaxy parameters
basic_map:
  distance: 1000000.0           # Distance in pc (e.g., 1e6 = 1 Mpc)
  inclination: 0.0              # Inclination in degrees
  posangle: 0.0                 # Position angle in degrees
  centrex: 512                  # Galaxy center X pixel (0-indexed)
  centrey: 512                  # Galaxy center Y pixel (0-indexed)
  minradius: 0.0                # Minimum radius in pc
  maxradius: 10000.0            # Maximum radius in pc

# Aperture settings
aperture:
  lapmin: 25.0                  # Minimum aperture diameter in pc
  lapmax: 6400.0                # Maximum aperture diameter in pc
  naperture: 9                  # Number of aperture sizes
  peak_res: 1                   # Index for peak identification
  max_res: 8                    # Index for maximum scale

# Peak identification
peak_id:
  npixmin: 20                   # Minimum pixels per peak
  nsigma: 5.0                   # Significance threshold
  logrange_s: 2.0               # Star contour log range
  logspacing_s: 0.5             # Star contour log spacing
  logrange_g: 2.0               # Gas contour log range
  logspacing_g: 0.5             # Gas contour log spacing

# Timeline parameters
timeline:
  tstariso: 1.0                 # Reference timescale in Myr
  tstariso_errmin: 0.0          # Downward error in Myr
  tstariso_errmax: 0.0          # Upward error in Myr
  tgasmini: 0.1                 # Minimum tgas in Myr
  tgasmaxi: 1000.0              # Maximum tgas in Myr
  tovermini: 0.01               # Minimum tover in Myr

# Fitting parameters
fitting:
  nmc: 1000                     # Monte Carlo samples
  ndepth: 4                     # Refinement depth
  ntry: 101                     # Grid size
  nphysmc: 1000000              # MC for derived physics

# Conversion factors (log10 values)
conversion:
  convstar: -3.69               # Star map -> SFR [Msun/yr]
  convgas: 2.31                 # Gas map -> mass [Msun]

# Analysis flags (true/false)
flags1:
  mask_images: true
  regrid: true
  smoothen: true
  sensitivity: true
  id_peaks: 1                   # 0=reuse, 1=identify, 2=reuse from path
  calc_ap_flux: true
  generate_plot: true
  get_distances: true
  calc_obs: true
  calc_fit: true
  diffuse_frac: true
  derive_phys: true
  write_output: true
  cleanup: 2                    # 0=keep, 1=prompt, 2=auto-delete
  autoexit: false

flags4:
  tophat: true                  # Use tophat (true) or Gaussian (false) kernel
  loglevels: true               # Log (true) or linear (false) contour spacing
  peak_prof: 2                  # 0=points, 1=discs, 2=Gaussians
  map_units: 1                  # 0=unknown, 1=SFR/gas, 2=gas/gas, 3=SFR/SFR
'''


def write_yaml_template(filepath: str | Path) -> None:
    """
    Write a YAML configuration template to a file.

    Args:
        filepath: Path where the template will be written
    """
    filepath = Path(filepath)
    filepath.write_text(YAML_TEMPLATE)

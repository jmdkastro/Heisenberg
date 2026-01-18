"""
Heisenberg configuration parameters.

This module defines the HeisenbergConfig class containing all ~180 parameters
needed for the Heisenberg analysis. Parameters are organized into logical groups
matching the IDL input_file structure.

The configuration can be loaded from an IDL-format input file or created
programmatically.
"""

from pathlib import Path
from typing import Optional, Literal
from pydantic import BaseModel, Field, field_validator, model_validator


class FileNames(BaseModel):
    """File name parameters for input data."""

    datadir: Path = Field(
        ...,
        description="Full path of data directory containing FITS files"
    )
    galaxy: str = Field(
        ...,
        description="Name of data set (used in output filenames)"
    )
    starfile: str = Field(
        ...,
        description="Name of primary stellar/SF map FITS file"
    )
    starfile2: Optional[str] = Field(
        None,
        description="Secondary stellar map for peak identification (if use_star2=1)"
    )
    gasfile: str = Field(
        ...,
        description="Name of primary gas map FITS file"
    )
    gasfile2: Optional[str] = Field(
        None,
        description="Secondary gas map for peak identification (if use_gas2=1)"
    )
    starfile3: Optional[str] = Field(
        None,
        description="Additional stellar map for flux ratio masking (if use_star3=1)"
    )


class PeakIdFileNames(BaseModel):
    """Peak identification file name parameters."""

    peaksdir: Optional[Path] = Field(
        None,
        description="Directory containing .sav files for identified peaks (if id_peaks=2)"
    )
    peakiddir: Optional[Path] = Field(
        None,
        description="Directory containing clumpfind output files (if id_peaks=2)"
    )
    starpeakidfile: Optional[str] = Field(
        None,
        description="Peak ID .dat file for starfile (if id_peaks=2)"
    )
    gaspeakidfile: Optional[str] = Field(
        None,
        description="Peak ID .dat file for gasfile (if id_peaks=2)"
    )
    intpeakidfile: Optional[str] = Field(
        None,
        description="Interactive peak ID report file (if id_peaks=2)"
    )


class MaskFileNames(BaseModel):
    """Mask file name parameters."""

    maskdir: Optional[Path] = Field(
        None,
        description="Full path of mask directory containing DS9 region files"
    )
    star_ext_mask: Optional[str] = Field(
        None,
        description="DS9 region file for external masking of starfile"
    )
    star_int_mask: Optional[str] = Field(
        None,
        description="DS9 region file for internal masking of starfile"
    )
    gas_ext_mask: Optional[str] = Field(
        None,
        description="DS9 region file for external masking of gasfile"
    )
    gas_int_mask: Optional[str] = Field(
        None,
        description="DS9 region file for internal masking of gasfile"
    )
    star_ext_mask2: Optional[str] = Field(
        None,
        description="DS9 region file for external masking of starfile2"
    )
    star_int_mask2: Optional[str] = Field(
        None,
        description="DS9 region file for internal masking of starfile2"
    )
    gas_ext_mask2: Optional[str] = Field(
        None,
        description="DS9 region file for external masking of gasfile2"
    )
    gas_int_mask2: Optional[str] = Field(
        None,
        description="DS9 region file for internal masking of gasfile2"
    )
    star_ext_mask3: Optional[str] = Field(
        None,
        description="DS9 region file for external masking of starfile3"
    )
    star_int_mask3: Optional[str] = Field(
        None,
        description="DS9 region file for internal masking of starfile3"
    )


class UnfilteredFileNames(BaseModel):
    """Unfiltered file name parameters for diffuse fraction calculation."""

    unfiltdir: Optional[Path] = Field(
        None,
        description="Directory containing unfiltered images"
    )
    star_unfilt_file: Optional[str] = Field(
        None,
        description="Unfiltered stellar map for diffuse fraction calculation"
    )
    gas_unfilt_file: Optional[str] = Field(
        None,
        description="Unfiltered gas map for diffuse fraction calculation"
    )


class Flags1(BaseModel):
    """Module switch flags (FLAGS 1)."""

    mask_images: bool = Field(
        True,
        description="Apply masks to images"
    )
    regrid: bool = Field(
        True,
        description="Read and regrid original files, synchronize masks across all images"
    )
    smoothen: bool = Field(
        True,
        description="Create smoothed maps for each aperture size"
    )
    sensitivity: bool = Field(
        True,
        description="Fit Gaussians to pixel intensity PDFs for sensitivity limits"
    )
    id_peaks: Literal[0, 1, 2] = Field(
        1,
        description="Peak identification mode: 0=reuse from default location, 1=identify new, 2=reuse from specified location"
    )
    calc_ap_flux: bool = Field(
        True,
        description="Calculate enclosed flux for each peak and aperture size"
    )
    generate_plot: bool = Field(
        True,
        description="Generate output plots and save to PostScript files"
    )
    get_distances: bool = Field(
        True,
        description="Calculate distances between all peak pairs"
    )
    calc_obs: bool = Field(
        True,
        description="Monte-Carlo sample peaks, get observed tuning fork"
    )
    calc_fit: bool = Field(
        True,
        description="Fit KL14 principle model"
    )
    diffuse_frac: bool = Field(
        True,
        description="Calculate diffuse fraction in images"
    )
    derive_phys: bool = Field(
        True,
        description="Calculate derived physical quantities"
    )
    write_output: bool = Field(
        True,
        description="Write results to output files"
    )
    cleanup: Literal[0, 1, 2] = Field(
        2,
        description="Cleanup mode: 0=keep files, 1=prompt before delete, 2=auto-delete"
    )
    autoexit: bool = Field(
        False,
        description="Exit automatically upon successful completion"
    )


class Flags2(BaseModel):
    """Ancillary file usage flags (FLAGS 2)."""

    use_star2: bool = Field(
        False,
        description="Use different map for identifying stellar peaks"
    )
    use_gas2: bool = Field(
        False,
        description="Use different map for identifying gas peaks"
    )
    use_star3: bool = Field(
        False,
        description="Use additional stellar tracer map for flux ratio masking"
    )


class Flags3(BaseModel):
    """Masking-related flags (FLAGS 3)."""

    mstar_ext: bool = Field(
        True,
        description="Mask starfile exterior to star_ext_mask regions"
    )
    mstar_int: bool = Field(
        True,
        description="Mask starfile interior to star_int_mask regions"
    )
    mgas_ext: bool = Field(
        True,
        description="Mask gasfile exterior to gas_ext_mask regions"
    )
    mgas_int: bool = Field(
        True,
        description="Mask gasfile interior to gas_int_mask regions"
    )
    mstar_ext2: bool = Field(
        False,
        description="Mask starfile2 exterior to star_ext_mask2 regions"
    )
    mstar_int2: bool = Field(
        False,
        description="Mask starfile2 interior to star_int_mask2 regions"
    )
    mgas_ext2: bool = Field(
        False,
        description="Mask gasfile2 exterior to gas_ext_mask2 regions"
    )
    mgas_int2: bool = Field(
        False,
        description="Mask gasfile2 interior to gas_int_mask2 regions"
    )
    mstar_ext3: bool = Field(
        False,
        description="Mask starfile3 exterior to star_ext_mask3 regions"
    )
    mstar_int3: bool = Field(
        False,
        description="Mask starfile3 interior to star_int_mask3 regions"
    )
    convert_masks: bool = Field(
        False,
        description="Convert masks from other DS9 coordinate systems to image coordinates"
    )
    cut_radius: bool = Field(
        True,
        description="Mask maps outside specified radial interval within galaxy"
    )


class Flags4(BaseModel):
    """Analysis option flags (FLAGS 4)."""

    set_centre: bool = Field(
        True,
        description="If True, use specified centrex/centrey; if False, use image center"
    )
    tophat: bool = Field(
        True,
        description="Use tophat kernel (True) or Gaussian kernel (False) for smoothing"
    )
    loglevels: bool = Field(
        True,
        description="Use logarithmic (True) or linear (False) contour level spacing"
    )
    peak_find_tui: bool = Field(
        False,
        description="Enable interactive text user interface for peak finding"
    )
    flux_weight: bool = Field(
        False,
        description="Use flux-weighted mean position (True) or brightest pixel (False)"
    )
    calc_ap_area: bool = Field(
        True,
        description="Calculate aperture area based on unmasked pixels"
    )
    tstar_incl: bool = Field(
        False,
        description="If True, tstariso includes overlap phase; if False, tstar=tstariso+tover"
    )
    peak_prof: Literal[0, 1, 2] = Field(
        2,
        description="Peak profile model: 0=points, 1=constant-density discs, 2=2D Gaussians"
    )
    map_units: Literal[0, 1, 2, 3] = Field(
        1,
        description="Map units: 0=unknown, 1=star=SFR/gas=gas, 2=star=gas1/gas=gas2, 3=star=SFR1/gas=SFR2"
    )
    star_tot_mode: Literal[0, 1] = Field(
        0,
        description="Star total mode: 0=calculate from map, 1=use star_tot_val"
    )
    gas_tot_mode: Literal[0, 1] = Field(
        0,
        description="Gas total mode: 0=calculate from map, 1=use gas_tot_val"
    )
    use_X11: bool = Field(
        True,
        description="Allow creation of X11 windows"
    )
    log10_output: bool = Field(
        True,
        description="Write output values as log10(value)"
    )


class BasicMapParams(BaseModel):
    """Basic map analysis parameters (INPUT PARAMETERS 1)."""

    distance: float = Field(
        ...,
        gt=0,
        description="Distance to galaxy in pc"
    )
    inclination: float = Field(
        0.0,
        ge=0,
        le=90,
        description="Inclination angle in degrees"
    )
    posangle: float = Field(
        0.0,
        description="Position angle in degrees"
    )
    centrex: int = Field(
        ...,
        ge=0,
        description="X-axis pixel coordinate of galaxy centre (0-indexed from left)"
    )
    centrey: int = Field(
        ...,
        ge=0,
        description="Y-axis pixel coordinate of galaxy centre (0-indexed from bottom)"
    )
    minradius: float = Field(
        0.0,
        ge=0,
        description="Minimum radius for analysis in pc"
    )
    maxradius: float = Field(
        ...,
        gt=0,
        description="Maximum radius for analysis in pc"
    )
    Fs1_Fs2_min: float = Field(
        15.0,
        description="Minimum primary-to-secondary SF tracer flux ratio (if use_star3=1)"
    )
    max_sample: int = Field(
        10,
        gt=0,
        description="Maximum number of pixels per map resolution FWHM"
    )
    astr_tolerance: float = Field(
        1.0e-6,
        gt=0,
        description="Allowable astrometric tolerance in decimal degrees"
    )
    nbins: int = Field(
        20,
        gt=0,
        description="Number of bins for sensitivity limit PDF fitting"
    )


class ApertureParams(BaseModel):
    """Aperture parameters (INPUT PARAMETERS 2)."""

    lapmin: float = Field(
        25.0,
        gt=0,
        description="Minimum aperture diameter in pc"
    )
    lapmax: float = Field(
        6400.0,
        gt=0,
        description="Maximum aperture diameter in pc"
    )
    naperture: int = Field(
        9,
        gt=0,
        description="Number of aperture sizes"
    )
    peak_res: int = Field(
        1,
        ge=0,
        description="Index of aperture size for peak identification and minimum fitting (0-indexed)"
    )
    max_res: int = Field(
        8,
        ge=0,
        description="Index of aperture size for maximum fitting and galactic averages (0-indexed)"
    )

    @model_validator(mode='after')
    def validate_aperture_indices(self) -> 'ApertureParams':
        if self.peak_res >= self.naperture:
            raise ValueError(f"peak_res ({self.peak_res}) must be < naperture ({self.naperture})")
        if self.max_res >= self.naperture:
            raise ValueError(f"max_res ({self.max_res}) must be < naperture ({self.naperture})")
        if self.peak_res > self.max_res:
            raise ValueError(f"peak_res ({self.peak_res}) must be <= max_res ({self.max_res})")
        return self


class PeakIdParams(BaseModel):
    """Peak identification parameters (INPUT PARAMETERS 3)."""

    npixmin: int = Field(
        20,
        ge=1,
        description="Minimum number of pixels for a valid peak"
    )
    nsigma: float = Field(
        5.0,
        gt=0,
        description="Multiple of sensitivity limit for valid peak"
    )
    logrange_s: float = Field(
        2.0,
        gt=0,
        description="Logarithmic range for stellar peak contour levels"
    )
    logspacing_s: float = Field(
        0.5,
        gt=0,
        description="Logarithmic interval between stellar contour levels"
    )
    logrange_g: float = Field(
        2.0,
        gt=0,
        description="Logarithmic range for gas peak contour levels"
    )
    logspacing_g: float = Field(
        0.5,
        gt=0,
        description="Logarithmic interval between gas contour levels"
    )
    nlinlevel_s: int = Field(
        11,
        gt=0,
        description="Number of linear contour levels for stellar peaks (if loglevels=0)"
    )
    nlinlevel_g: int = Field(
        11,
        gt=0,
        description="Number of linear contour levels for gas peaks (if loglevels=0)"
    )


class TimelineParams(BaseModel):
    """Timeline parameters (INPUT PARAMETERS 4)."""

    tstariso: float = Field(
        1.0,
        gt=0,
        description="Reference timescale of star formation tracer in Myr"
    )
    tstariso_errmin: float = Field(
        0.0,
        ge=0,
        description="Downward standard error of tstariso in Myr"
    )
    tstariso_errmax: float = Field(
        0.0,
        ge=0,
        description="Upward standard error of tstariso in Myr"
    )
    tgasmini: float = Field(
        0.1,
        gt=0,
        description="Minimum tgas value during fitting in Myr"
    )
    tgasmaxi: float = Field(
        1000.0,
        gt=0,
        description="Maximum tgas value during fitting in Myr"
    )
    tovermini: float = Field(
        0.01,
        gt=0,
        description="Minimum tover value during fitting in Myr"
    )


class FittingParams(BaseModel):
    """Fitting parameters (INPUT PARAMETERS 5)."""

    nmc: int = Field(
        1000,
        gt=0,
        description="Number of Monte Carlo peak drawing experiments"
    )
    ndepth: int = Field(
        4,
        gt=0,
        description="Maximum number of parameter refinement loops"
    )
    ntry: int = Field(
        101,
        gt=0,
        description="Size of each parameter array for grid search"
    )
    nphysmc: int = Field(
        1000000,
        gt=0,
        description="Number of Monte Carlo experiments for derived physics error propagation"
    )


class FourierFilterParams(BaseModel):
    """Fourier filtering parameters (INPUT PARAMETERS 6)."""

    use_unfilt_ims: bool = Field(
        False,
        description="Calculate diffuse fraction using unfiltered images"
    )
    diffuse_quant: Literal[0, 1] = Field(
        1,
        description="Diffuse quantity: 0=flux, 1=power"
    )
    f_filter_type: Literal[0, 1, 2] = Field(
        2,
        description="Filter type: 0=butterworth, 1=gaussian, 2=ideal"
    )
    bw_order: int = Field(
        2,
        ge=1,
        description="Butterworth filter order"
    )
    filter_len_conv: float = Field(
        1.0,
        gt=0,
        description="Conversion factor for Fourier filter cut length (cut_length = lambda * filter_len_conv)"
    )
    emfrac_cor_mode: Literal[0, 1, 2, 3, 4, 5] = Field(
        0,
        description="Emission fraction correction mode: 0=none, 1=flux loss, 4=overlap, 5=both"
    )
    rpeak_cor_mode: Literal[0, 1] = Field(
        0,
        description="rpeak correction mode: 0=measured, 1=supplied values"
    )
    rpeaks_cor_val: float = Field(
        1.0,
        gt=0,
        description="r_peak_star value for flux loss correction (if rpeak_cor_mode=1)"
    )
    rpeaks_cor_emin: float = Field(
        0.5,
        ge=0,
        description="Downward error of r_peak_star"
    )
    rpeaks_cor_emax: float = Field(
        0.5,
        ge=0,
        description="Upward error of r_peak_star"
    )
    rpeakg_cor_val: float = Field(
        1.0,
        gt=0,
        description="r_peak_gas value for flux loss correction (if rpeak_cor_mode=1)"
    )
    rpeakg_cor_emin: float = Field(
        0.5,
        ge=0,
        description="Downward error of r_peak_gas"
    )
    rpeakg_cor_emax: float = Field(
        0.5,
        ge=0,
        description="Upward error of r_peak_gas"
    )


class ConversionParams(BaseModel):
    """Conversion and constant parameters (INPUT PARAMETERS 7)."""

    convstar: float = Field(
        -3.69206,
        description="Log of conversion factor: pixel value -> SFR (Msun/yr) or gas mass (Msun)"
    )
    convstar_rerr: float = Field(
        0.0,
        ge=0,
        description="Relative error of convstar"
    )
    convgas: float = Field(
        2.30794,
        description="Log of conversion factor: pixel value -> gas mass (Msun) or SFR (Msun/yr)"
    )
    convgas_rerr: float = Field(
        0.0,
        ge=0,
        description="Relative error of convgas"
    )
    convstar3: float = Field(
        0.0,
        description="Log of conversion factor for starfile3"
    )
    convstar3_rerr: float = Field(
        0.0,
        ge=0,
        description="Relative error of convstar3"
    )
    lighttomass: float = Field(
        0.002,
        description="Light-to-mass ratio of feedback mechanism in m^2 s^-3 (default: SNe)"
    )
    momratetomass: float = Field(
        5.0e-10,
        description="Momentum output rate per unit mass in m s^-2 (default: SNe+winds)"
    )
    star_tot_val: float = Field(
        1.0,
        description="Input value of total SFR/gas mass for stellar map (if star_tot_mode=1)"
    )
    star_tot_err: float = Field(
        0.1,
        ge=0,
        description="Error on star_tot_val"
    )
    gas_tot_val: float = Field(
        1.0e9,
        description="Input value of total gas mass/SFR for gas map (if gas_tot_mode=1)"
    )
    gas_tot_err: float = Field(
        1.0e8,
        ge=0,
        description="Error on gas_tot_val"
    )


class SensitivityParams(BaseModel):
    """Sensitivity parameters (INPUT PARAMETERS 8)."""

    use_stds: bool = Field(
        False,
        description="Use supplied standard deviations (True) or calculate (False)"
    )
    std_star: float = Field(
        0.1,
        gt=0,
        description="Standard deviation of starfile (if use_stds=True)"
    )
    std_star3: float = Field(
        0.1,
        gt=0,
        description="Standard deviation of starfile3 (if use_stds=True)"
    )
    std_gas: float = Field(
        0.1,
        gt=0,
        description="Standard deviation of gasfile (if use_stds=True)"
    )


class NoiseThresholdParams(BaseModel):
    """Noise threshold parameters (INPUT PARAMETERS 9)."""

    use_noisecut: bool = Field(
        True,
        description="Mask values below noise threshold after filtering"
    )
    noisethresh_s: float = Field(
        20.0,
        description="Noise threshold for star map after filtering"
    )
    noisethresh_g: float = Field(
        20.0,
        description="Noise threshold for gas map after filtering"
    )


class IterationParams(BaseModel):
    """Fourier diffuse removal iteration parameters (INPUT PARAMETERS 10)."""

    use_guess: bool = Field(
        False,
        description="Filter images with initial lambda guess before first run"
    )
    initial_guess: float = Field(
        200.0,
        gt=0,
        description="Initial estimate of lambda in pc for pre-filtering"
    )
    iter_criterion: float = Field(
        0.05,
        gt=0,
        lt=1,
        description="Fractional convergence criterion for lambda"
    )
    iter_crit_len: int = Field(
        2,
        ge=1,
        description="Number of previous iterations for convergence check"
    )
    iter_nmax: int = Field(
        10,
        ge=1,
        description="Maximum number of iterations"
    )
    iter_filter: Literal[0, 1, 2] = Field(
        2,
        description="Iteration filter type: 0=butterworth, 1=gaussian, 2=ideal"
    )
    iter_bwo: int = Field(
        2,
        ge=1,
        description="Butterworth order for iteration filtering"
    )
    iter_len_conv: float = Field(
        2.0,
        gt=0,
        description="Conversion factor for iteration filter cut length"
    )
    iter_rpeak_mode: Literal[0, 1] = Field(
        0,
        description="Iteration rpeak mode: 0=use rpeak_cor_mode, 1=use iter0 values"
    )
    iter_tot_mode_s: Literal[0, 1, 2] = Field(
        0,
        description="Iteration star total mode: 0=per-iteration, 1=from iter0, 2=from star_tot_val"
    )
    iter_tot_mode_g: Literal[0, 1, 2] = Field(
        0,
        description="Iteration gas total mode: 0=per-iteration, 1=from iter0, 2=from gas_tot_val"
    )
    iter_autoexit: bool = Field(
        False,
        description="Exit automatically upon successful iteration completion"
    )
    use_nice: bool = Field(
        False,
        description="Use nice command for spawned processes"
    )
    nice_value: int = Field(
        0,
        ge=-20,
        le=19,
        description="Nice priority value (-20=highest, 19=lowest)"
    )


class HeisenbergConfig(BaseModel):
    """
    Complete configuration for a Heisenberg analysis run.

    This class contains all ~180 parameters organized into logical groups.
    It can be loaded from an IDL-format input file or constructed programmatically.

    Example:
        >>> from heisenberg.config import HeisenbergConfig, load_config
        >>> config = load_config("input_file")
        >>> print(config.files.galaxy)
        >>> print(config.basic_map.distance)
    """

    # File names
    files: FileNames
    peak_files: PeakIdFileNames = Field(default_factory=PeakIdFileNames)
    mask_files: MaskFileNames = Field(default_factory=MaskFileNames)
    unfilt_files: UnfilteredFileNames = Field(default_factory=UnfilteredFileNames)

    # Flags
    flags1: Flags1 = Field(default_factory=Flags1)
    flags2: Flags2 = Field(default_factory=Flags2)
    flags3: Flags3 = Field(default_factory=Flags3)
    flags4: Flags4 = Field(default_factory=Flags4)

    # Input parameters
    basic_map: BasicMapParams
    aperture: ApertureParams = Field(default_factory=ApertureParams)
    peak_id: PeakIdParams = Field(default_factory=PeakIdParams)
    timeline: TimelineParams = Field(default_factory=TimelineParams)
    fitting: FittingParams = Field(default_factory=FittingParams)
    fourier: FourierFilterParams = Field(default_factory=FourierFilterParams)
    conversion: ConversionParams = Field(default_factory=ConversionParams)
    sensitivity: SensitivityParams = Field(default_factory=SensitivityParams)
    noise: NoiseThresholdParams = Field(default_factory=NoiseThresholdParams)
    iteration: IterationParams = Field(default_factory=IterationParams)

    model_config = {
        "extra": "forbid",  # Don't allow extra fields
        "validate_assignment": True,  # Validate on attribute assignment
    }

    @classmethod
    def from_file(cls, filepath: str | Path) -> "HeisenbergConfig":
        """
        Load configuration from an IDL-format input file.

        Args:
            filepath: Path to the input file

        Returns:
            HeisenbergConfig instance
        """
        from heisenberg.config.reader import load_config
        return load_config(filepath)

    # Convenience properties for flat access to commonly used parameters
    @property
    def galaxy(self) -> str:
        """Galaxy name for output files."""
        return self.files.galaxy

    @property
    def distance(self) -> float:
        """Distance to galaxy in pc."""
        return self.basic_map.distance

    @property
    def datadir(self) -> Path:
        """Data directory path."""
        return self.files.datadir

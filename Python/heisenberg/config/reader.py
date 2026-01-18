"""
Input file parser for Heisenberg configuration.

This module provides functions to parse configuration files into
HeisenbergConfig instances. Supports two formats:

1. IDL format (.inp, no extension, or any non-.yaml/.yml file):
   - Space-separated key-value pairs
   - # comments (line and inline)
   - "-" as null/None marker for optional paths
   - Automatic type conversion (int, float, string)

2. YAML format (.yaml, .yml):
   - Standard YAML syntax
   - Nested structure matching HeisenbergConfig layout
   - null for optional paths
"""

import re
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

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


def parse_value(value_str: str) -> Any:
    """
    Parse a string value into the appropriate Python type.

    Args:
        value_str: String value from input file

    Returns:
        Parsed value (int, float, bool, str, or None)
    """
    # Handle null marker
    if value_str == "-":
        return None

    # Try integer
    try:
        # Check if it looks like an integer (no decimal point or exponent)
        if re.match(r'^-?\d+$', value_str):
            return int(value_str)
    except ValueError:
        pass

    # Try float
    try:
        return float(value_str)
    except ValueError:
        pass

    # Return as string
    return value_str


def parse_input_file(filepath: str | Path) -> Dict[str, Any]:
    """
    Parse an IDL-format input file into a dictionary.

    The format is:
    - One parameter per line
    - Parameter name followed by value, separated by whitespace
    - Lines starting with # are comments
    - Inline comments after # are stripped
    - "-" represents null/disabled value

    Args:
        filepath: Path to the input file

    Returns:
        Dictionary of parameter name -> value
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    params: Dict[str, Any] = {}

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Strip whitespace
            line = line.strip()

            # Skip empty lines and comment lines
            if not line or line.startswith('#'):
                continue

            # Remove inline comments
            if '#' in line:
                line = line[:line.index('#')].strip()

            # Skip if nothing left after removing comments
            if not line:
                continue

            # Split on whitespace
            parts = line.split()
            if len(parts) < 2:
                continue  # Need at least name and value

            param_name = parts[0]
            param_value = parts[1]

            # Parse and store
            params[param_name] = parse_value(param_value)

    return params


def parse_yaml_file(filepath: str | Path) -> Dict[str, Any]:
    """
    Parse a YAML configuration file into a flat dictionary.

    The YAML file should have a nested structure matching the
    HeisenbergConfig layout. This function flattens it into
    a dictionary compatible with the IDL-format parameter names.

    Args:
        filepath: Path to the YAML file

    Returns:
        Dictionary of parameter name -> value (flattened)
    """
    if not YAML_AVAILABLE:
        raise ImportError(
            "PyYAML is required to parse YAML configuration files. "
            "Install it with: pip install pyyaml"
        )

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Configuration file not found: {filepath}")

    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

    if data is None:
        return {}

    # Flatten nested structure to match IDL parameter names
    params: Dict[str, Any] = {}

    # Map YAML section names to their parameters
    section_mapping = {
        'files': [
            'datadir', 'galaxy', 'starfile', 'starfile2', 'gasfile',
            'gasfile2', 'starfile3'
        ],
        'peak_files': [
            'peaksdir', 'peakiddir', 'starpeakidfile', 'gaspeakidfile',
            'intpeakidfile'
        ],
        'mask_files': [
            'maskdir', 'star_ext_mask', 'star_int_mask', 'gas_ext_mask',
            'gas_int_mask', 'star_ext_mask2', 'star_int_mask2',
            'gas_ext_mask2', 'gas_int_mask2', 'star_ext_mask3', 'star_int_mask3'
        ],
        'unfilt_files': [
            'unfiltdir', 'star_unfilt_file', 'gas_unfilt_file'
        ],
        'flags1': [
            'mask_images', 'regrid', 'smoothen', 'sensitivity', 'id_peaks',
            'calc_ap_flux', 'generate_plot', 'get_distances', 'calc_obs',
            'calc_fit', 'diffuse_frac', 'derive_phys', 'write_output',
            'cleanup', 'autoexit'
        ],
        'flags2': ['use_star2', 'use_gas2', 'use_star3'],
        'flags3': [
            'mstar_ext', 'mstar_int', 'mgas_ext', 'mgas_int',
            'mstar_ext2', 'mstar_int2', 'mgas_ext2', 'mgas_int2',
            'mstar_ext3', 'mstar_int3', 'convert_masks', 'cut_radius'
        ],
        'flags4': [
            'set_centre', 'tophat', 'loglevels', 'peak_find_tui',
            'flux_weight', 'calc_ap_area', 'tstar_incl', 'peak_prof',
            'map_units', 'star_tot_mode', 'gas_tot_mode', 'use_X11',
            'log10_output'
        ],
        'basic_map': [
            'distance', 'inclination', 'posangle', 'centrex', 'centrey',
            'minradius', 'maxradius', 'Fs1_Fs2_min', 'max_sample',
            'astr_tolerance', 'nbins'
        ],
        'aperture': [
            'lapmin', 'lapmax', 'naperture', 'peak_res', 'max_res'
        ],
        'peak_id': [
            'npixmin', 'nsigma', 'logrange_s', 'logspacing_s',
            'logrange_g', 'logspacing_g', 'nlinlevel_s', 'nlinlevel_g'
        ],
        'timeline': [
            'tstariso', 'tstariso_errmin', 'tstariso_errmax',
            'tgasmini', 'tgasmaxi', 'tovermini'
        ],
        'fitting': ['nmc', 'ndepth', 'ntry', 'nphysmc'],
        'fourier': [
            'use_unfilt_ims', 'diffuse_quant', 'f_filter_type', 'bw_order',
            'filter_len_conv', 'emfrac_cor_mode', 'rpeak_cor_mode',
            'rpeaks_cor_val', 'rpeaks_cor_emin', 'rpeaks_cor_emax',
            'rpeakg_cor_val', 'rpeakg_cor_emin', 'rpeakg_cor_emax'
        ],
        'conversion': [
            'convstar', 'convstar_rerr', 'convgas', 'convgas_rerr',
            'convstar3', 'convstar3_rerr', 'lighttomass', 'momratetomass',
            'star_tot_val', 'star_tot_err', 'gas_tot_val', 'gas_tot_err'
        ],
        'sensitivity': ['use_stds', 'std_star', 'std_star3', 'std_gas'],
        'noise': ['use_noisecut', 'noisethresh_s', 'noisethresh_g'],
        'iteration': [
            'use_guess', 'initial_guess', 'iter_criterion', 'iter_crit_len',
            'iter_nmax', 'iter_filter', 'iter_bwo', 'iter_len_conv',
            'iter_rpeak_mode', 'iter_tot_mode_s', 'iter_tot_mode_g',
            'iter_autoexit', 'use_nice', 'nice_value'
        ],
    }

    # Extract parameters from each section
    for section, keys in section_mapping.items():
        section_data = data.get(section, {})
        if section_data is None:
            section_data = {}
        for key in keys:
            if key in section_data:
                value = section_data[key]
                # Convert YAML booleans to integers for flags (IDL compatibility)
                if isinstance(value, bool) and section.startswith('flags'):
                    value = 1 if value else 0
                params[key] = value

    return params


def _get_optional_path(params: Dict[str, Any], key: str) -> Optional[Path]:
    """Get an optional path parameter, converting string to Path."""
    value = params.get(key)
    if value is None or value == "-":
        return None
    return Path(value)


def _get_optional_str(params: Dict[str, Any], key: str) -> Optional[str]:
    """Get an optional string parameter."""
    value = params.get(key)
    if value is None or value == "-":
        return None
    return str(value)


def _get_bool(params: Dict[str, Any], key: str, default: bool = False) -> bool:
    """Get a boolean parameter (0/1 -> False/True)."""
    value = params.get(key, default)
    if isinstance(value, bool):
        return value
    return bool(int(value))


def _get_int(params: Dict[str, Any], key: str, default: int = 0) -> int:
    """Get an integer parameter."""
    value = params.get(key, default)
    return int(value)


def _get_float(params: Dict[str, Any], key: str, default: float = 0.0) -> float:
    """Get a float parameter."""
    value = params.get(key, default)
    return float(value)


def _is_yaml_file(filepath: Path) -> bool:
    """Check if a file is a YAML file based on extension."""
    return filepath.suffix.lower() in ('.yaml', '.yml')


def load_config(filepath: str | Path) -> HeisenbergConfig:
    """
    Load a HeisenbergConfig from an input file.

    Supports two formats:
    - YAML (.yaml, .yml): Modern format with nested structure
    - IDL format (any other extension): Legacy format from IDL codebase

    Args:
        filepath: Path to the configuration file

    Returns:
        HeisenbergConfig instance with all parameters populated

    Raises:
        FileNotFoundError: If the configuration file doesn't exist
        ImportError: If PyYAML is not installed for YAML files
    """
    filepath = Path(filepath)

    if _is_yaml_file(filepath):
        params = parse_yaml_file(filepath)
    else:
        params = parse_input_file(filepath)

    # Build file names
    files = FileNames(
        datadir=Path(params['datadir']),
        galaxy=str(params['galaxy']),
        starfile=str(params['starfile']),
        starfile2=_get_optional_str(params, 'starfile2'),
        gasfile=str(params['gasfile']),
        gasfile2=_get_optional_str(params, 'gasfile2'),
        starfile3=_get_optional_str(params, 'starfile3'),
    )

    # Build peak ID file names
    peak_files = PeakIdFileNames(
        peaksdir=_get_optional_path(params, 'peaksdir'),
        peakiddir=_get_optional_path(params, 'peakiddir'),
        starpeakidfile=_get_optional_str(params, 'starpeakidfile'),
        gaspeakidfile=_get_optional_str(params, 'gaspeakidfile'),
        intpeakidfile=_get_optional_str(params, 'intpeakidfile'),
    )

    # Build mask file names
    mask_files = MaskFileNames(
        maskdir=_get_optional_path(params, 'maskdir'),
        star_ext_mask=_get_optional_str(params, 'star_ext_mask'),
        star_int_mask=_get_optional_str(params, 'star_int_mask'),
        gas_ext_mask=_get_optional_str(params, 'gas_ext_mask'),
        gas_int_mask=_get_optional_str(params, 'gas_int_mask'),
        star_ext_mask2=_get_optional_str(params, 'star_ext_mask2'),
        star_int_mask2=_get_optional_str(params, 'star_int_mask2'),
        gas_ext_mask2=_get_optional_str(params, 'gas_ext_mask2'),
        gas_int_mask2=_get_optional_str(params, 'gas_int_mask2'),
        star_ext_mask3=_get_optional_str(params, 'star_ext_mask3'),
        star_int_mask3=_get_optional_str(params, 'star_int_mask3'),
    )

    # Build unfiltered file names
    unfilt_files = UnfilteredFileNames(
        unfiltdir=_get_optional_path(params, 'unfiltdir'),
        star_unfilt_file=_get_optional_str(params, 'star_unfilt_file'),
        gas_unfilt_file=_get_optional_str(params, 'gas_unfilt_file'),
    )

    # Build flags
    flags1 = Flags1(
        mask_images=_get_bool(params, 'mask_images', True),
        regrid=_get_bool(params, 'regrid', True),
        smoothen=_get_bool(params, 'smoothen', True),
        sensitivity=_get_bool(params, 'sensitivity', True),
        id_peaks=_get_int(params, 'id_peaks', 1),
        calc_ap_flux=_get_bool(params, 'calc_ap_flux', True),
        generate_plot=_get_bool(params, 'generate_plot', True),
        get_distances=_get_bool(params, 'get_distances', True),
        calc_obs=_get_bool(params, 'calc_obs', True),
        calc_fit=_get_bool(params, 'calc_fit', True),
        diffuse_frac=_get_bool(params, 'diffuse_frac', True),
        derive_phys=_get_bool(params, 'derive_phys', True),
        write_output=_get_bool(params, 'write_output', True),
        cleanup=_get_int(params, 'cleanup', 2),
        autoexit=_get_bool(params, 'autoexit', False),
    )

    flags2 = Flags2(
        use_star2=_get_bool(params, 'use_star2', False),
        use_gas2=_get_bool(params, 'use_gas2', False),
        use_star3=_get_bool(params, 'use_star3', False),
    )

    flags3 = Flags3(
        mstar_ext=_get_bool(params, 'mstar_ext', True),
        mstar_int=_get_bool(params, 'mstar_int', True),
        mgas_ext=_get_bool(params, 'mgas_ext', True),
        mgas_int=_get_bool(params, 'mgas_int', True),
        mstar_ext2=_get_bool(params, 'mstar_ext2', False),
        mstar_int2=_get_bool(params, 'mstar_int2', False),
        mgas_ext2=_get_bool(params, 'mgas_ext2', False),
        mgas_int2=_get_bool(params, 'mgas_int2', False),
        mstar_ext3=_get_bool(params, 'mstar_ext3', False),
        mstar_int3=_get_bool(params, 'mstar_int3', False),
        convert_masks=_get_bool(params, 'convert_masks', False),
        cut_radius=_get_bool(params, 'cut_radius', True),
    )

    flags4 = Flags4(
        set_centre=_get_bool(params, 'set_centre', True),
        tophat=_get_bool(params, 'tophat', True),
        loglevels=_get_bool(params, 'loglevels', True),
        peak_find_tui=_get_bool(params, 'peak_find_tui', False),
        flux_weight=_get_bool(params, 'flux_weight', False),
        calc_ap_area=_get_bool(params, 'calc_ap_area', True),
        tstar_incl=_get_bool(params, 'tstar_incl', False),
        peak_prof=_get_int(params, 'peak_prof', 2),
        map_units=_get_int(params, 'map_units', 1),
        star_tot_mode=_get_int(params, 'star_tot_mode', 0),
        gas_tot_mode=_get_int(params, 'gas_tot_mode', 0),
        use_X11=_get_bool(params, 'use_X11', True),
        log10_output=_get_bool(params, 'log10_output', True),
    )

    # Build input parameters
    basic_map = BasicMapParams(
        distance=_get_float(params, 'distance'),
        inclination=_get_float(params, 'inclination', 0.0),
        posangle=_get_float(params, 'posangle', 0.0),
        centrex=_get_int(params, 'centrex'),
        centrey=_get_int(params, 'centrey'),
        minradius=_get_float(params, 'minradius', 0.0),
        maxradius=_get_float(params, 'maxradius'),
        Fs1_Fs2_min=_get_float(params, 'Fs1_Fs2_min', 15.0),
        max_sample=_get_int(params, 'max_sample', 10),
        astr_tolerance=_get_float(params, 'astr_tolerance', 1.0e-6),
        nbins=_get_int(params, 'nbins', 20),
    )

    aperture = ApertureParams(
        lapmin=_get_float(params, 'lapmin', 25.0),
        lapmax=_get_float(params, 'lapmax', 6400.0),
        naperture=_get_int(params, 'naperture', 9),
        peak_res=_get_int(params, 'peak_res', 1),
        max_res=_get_int(params, 'max_res', 8),
    )

    peak_id = PeakIdParams(
        npixmin=_get_int(params, 'npixmin', 20),
        nsigma=_get_float(params, 'nsigma', 5.0),
        logrange_s=_get_float(params, 'logrange_s', 2.0),
        logspacing_s=_get_float(params, 'logspacing_s', 0.5),
        logrange_g=_get_float(params, 'logrange_g', 2.0),
        logspacing_g=_get_float(params, 'logspacing_g', 0.5),
        nlinlevel_s=_get_int(params, 'nlinlevel_s', 11),
        nlinlevel_g=_get_int(params, 'nlinlevel_g', 11),
    )

    timeline = TimelineParams(
        tstariso=_get_float(params, 'tstariso', 1.0),
        tstariso_errmin=_get_float(params, 'tstariso_errmin', 0.0),
        tstariso_errmax=_get_float(params, 'tstariso_errmax', 0.0),
        tgasmini=_get_float(params, 'tgasmini', 0.1),
        tgasmaxi=_get_float(params, 'tgasmaxi', 1000.0),
        tovermini=_get_float(params, 'tovermini', 0.01),
    )

    fitting = FittingParams(
        nmc=_get_int(params, 'nmc', 1000),
        ndepth=_get_int(params, 'ndepth', 4),
        ntry=_get_int(params, 'ntry', 101),
        nphysmc=_get_int(params, 'nphysmc', 1000000),
    )

    fourier = FourierFilterParams(
        use_unfilt_ims=_get_bool(params, 'use_unfilt_ims', False),
        diffuse_quant=_get_int(params, 'diffuse_quant', 1),
        f_filter_type=_get_int(params, 'f_filter_type', 2),
        bw_order=_get_int(params, 'bw_order', 2),
        filter_len_conv=_get_float(params, 'filter_len_conv', 1.0),
        emfrac_cor_mode=_get_int(params, 'emfrac_cor_mode', 0),
        rpeak_cor_mode=_get_int(params, 'rpeak_cor_mode', 0),
        rpeaks_cor_val=_get_float(params, 'rpeaks_cor_val', 1.0),
        rpeaks_cor_emin=_get_float(params, 'rpeaks_cor_emin', 0.5),
        rpeaks_cor_emax=_get_float(params, 'rpeaks_cor_emax', 0.5),
        rpeakg_cor_val=_get_float(params, 'rpeakg_cor_val', 1.0),
        rpeakg_cor_emin=_get_float(params, 'rpeakg_cor_emin', 0.5),
        rpeakg_cor_emax=_get_float(params, 'rpeakg_cor_emax', 0.5),
    )

    conversion = ConversionParams(
        convstar=_get_float(params, 'convstar', -3.69206),
        convstar_rerr=_get_float(params, 'convstar_rerr', 0.0),
        convgas=_get_float(params, 'convgas', 2.30794),
        convgas_rerr=_get_float(params, 'convgas_rerr', 0.0),
        convstar3=_get_float(params, 'convstar3', 0.0),
        convstar3_rerr=_get_float(params, 'convstar3_rerr', 0.0),
        lighttomass=_get_float(params, 'lighttomass', 0.002),
        momratetomass=_get_float(params, 'momratetomass', 5.0e-10),
        star_tot_val=_get_float(params, 'star_tot_val', 1.0),
        star_tot_err=_get_float(params, 'star_tot_err', 0.1),
        gas_tot_val=_get_float(params, 'gas_tot_val', 1.0e9),
        gas_tot_err=_get_float(params, 'gas_tot_err', 1.0e8),
    )

    sensitivity = SensitivityParams(
        use_stds=_get_bool(params, 'use_stds', False),
        std_star=_get_float(params, 'std_star', 0.1),
        std_star3=_get_float(params, 'std_star3', 0.1),
        std_gas=_get_float(params, 'std_gas', 0.1),
    )

    noise = NoiseThresholdParams(
        use_noisecut=_get_bool(params, 'use_noisecut', True),
        noisethresh_s=_get_float(params, 'noisethresh_s', 20.0),
        noisethresh_g=_get_float(params, 'noisethresh_g', 20.0),
    )

    iteration = IterationParams(
        use_guess=_get_bool(params, 'use_guess', False),
        initial_guess=_get_float(params, 'initial_guess', 200.0),
        iter_criterion=_get_float(params, 'iter_criterion', 0.05),
        iter_crit_len=_get_int(params, 'iter_crit_len', 2),
        iter_nmax=_get_int(params, 'iter_nmax', 10),
        iter_filter=_get_int(params, 'iter_filter', 2),
        iter_bwo=_get_int(params, 'iter_bwo', 2),
        iter_len_conv=_get_float(params, 'iter_len_conv', 2.0),
        iter_rpeak_mode=_get_int(params, 'iter_rpeak_mode', 0),
        iter_tot_mode_s=_get_int(params, 'iter_tot_mode_s', 0),
        iter_tot_mode_g=_get_int(params, 'iter_tot_mode_g', 0),
        iter_autoexit=_get_bool(params, 'iter_autoexit', False),
        use_nice=_get_bool(params, 'use_nice', False),
        nice_value=_get_int(params, 'nice_value', 0),
    )

    # Build complete config
    return HeisenbergConfig(
        files=files,
        peak_files=peak_files,
        mask_files=mask_files,
        unfilt_files=unfilt_files,
        flags1=flags1,
        flags2=flags2,
        flags3=flags3,
        flags4=flags4,
        basic_map=basic_map,
        aperture=aperture,
        peak_id=peak_id,
        timeline=timeline,
        fitting=fitting,
        fourier=fourier,
        conversion=conversion,
        sensitivity=sensitivity,
        noise=noise,
        iteration=iteration,
    )

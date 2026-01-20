# Heisenberg

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

**Measuring the physics of star formation and feedback from observations of gas and young stars.**

## Overview

Heisenberg implements the "uncertainty principle for star formation" methodology developed by Kruijssen & Longmore (2014). By analysing the spatial decorrelation between gas and stellar tracers as a function of spatial scale, the code derives fundamental timescales governing the star formation process:

- **Gas cloud lifetime** (t_gas): The time from molecular cloud assembly to dispersal
- **Feedback timescale** (t_over): The duration of the embedded phase when both gas and young stars are visible
- **Region separation length** (lambda): The characteristic spacing between independent star-forming regions
- **Star formation efficiency**: The fraction of gas converted to stars

Both **Python** (recommended) and **IDL** (legacy) implementations are provided for full backward compatibility.

## Scientific Background

The method exploits the fact that gas and young stellar tracers are spatially correlated on large scales (where many independent regions are averaged) but decorrelated on small scales (where individual regions are resolved). The characteristic "tuning fork" shape of the flux ratio versus aperture size diagram encodes information about the underlying evolutionary timeline.

For detailed scientific background, see:
- [Scientific Background](docs/scientific-background.md)
- Kruijssen & Longmore (2014), MNRAS 439, 3239
- Kruijssen et al. (2018), MNRAS 479, 1866

## Key Features

- **Dual implementation**: Modern Python package with full test coverage + legacy IDL code
- **Two execution modes**:
  - Single-pass analysis (`heisenberg_nodf`)
  - Iterative diffuse filtering (`heisenberg`) for robust lambda determination
- **Interactive TUI**: Visual peak finding with real-time parameter tuning
- **Flexible I/O**: YAML or IDL-format configuration; JSON/text output; DS9 region export
- **Comprehensive plotting**: Tuning fork diagrams, parameter PDFs, iteration tracking

## Installation

### Python (Recommended)

```bash
# From source
cd Python/
pip install -e .

# Or with development dependencies
pip install -e ".[dev]"
```

**Requirements**: Python 3.9+, numpy, scipy, astropy, matplotlib, click

### IDL (Legacy)

**Requirements**: IDL 8.0+, [IDL Astronomy User's Library](https://idlastro.gsfc.nasa.gov/)

```idl
; Add to your IDL startup file or run manually:
!PATH = !PATH + ':' + expand_path('+/path/to/heisenberg/IDL/')
```

## Quick Start

### Python CLI

```bash
# Full analysis with iterative diffuse filtering (recommended)
heisenberg run config.yaml

# Single-pass without diffuse filtering
heisenberg run config.yaml --no-diffuse-filtering

# Interactive peak finding to tune detection parameters
heisenberg run config.yaml --interactive

# Validate configuration file
heisenberg validate config.yaml

# Quick test with synthetic data
heisenberg quick-test --plot
```

### Python API

```python
from heisenberg.core import run_tuningfork, TuningForkConfig
from heisenberg.io import read_fits

# Load data
star_data, star_wcs, _ = read_fits("star_map.fits")
gas_data, gas_wcs, _ = read_fits("gas_map.fits")

# Run analysis
config = TuningForkConfig(
    lap_min=25.0,      # Minimum aperture diameter (pc)
    lap_max=3200.0,    # Maximum aperture diameter (pc)
    n_apertures=9,     # Number of aperture sizes
    tstar=4.0,         # Stellar tracer timescale (Myr)
)

result = run_tuningfork(star_data, gas_data, pixel_scale=10.0, config=config)

print(f"t_gas  = {result.fit.tgas:.1f} (+{result.fit.tgas_errmax:.1f}/-{result.fit.tgas_errmin:.1f}) Myr")
print(f"t_over = {result.fit.tover:.1f} (+{result.fit.tover_errmax:.1f}/-{result.fit.tover_errmin:.1f}) Myr")
print(f"lambda = {result.fit.lambda_:.0f} (+{result.fit.lambda_errmax:.0f}/-{result.fit.lambda_errmin:.0f}) pc")
```

### IDL

```idl
; With iterative diffuse filtering (recommended)
idl heisenberg -arg /path/to/input_file

; Or interactively:
IDL> inputfile = '/path/to/input_file'
IDL> .run heisenberg
```

```idl
; Single-pass without diffuse filtering
idl heisenberg_nodf -arg /path/to/input_file
```

## Execution Modes

### 1. Iterative Diffuse Filtering (Default)

**Entry points**: `heisenberg run config.yaml` (Python) / `heisenberg.pro` (IDL)

Iteratively filters diffuse emission using the fitted lambda value, then refits until lambda converges. This is the recommended mode for most scientific analyses as it properly accounts for extended emission.

### 2. Single-Pass Analysis

**Entry points**: `heisenberg run config.yaml -n` (Python) / `heisenberg_nodf.pro` (IDL)

Runs a single tuningfork fit without diffuse filtering. Useful for:
- Quick tests and parameter exploration
- Pre-filtered data
- Comparing with/without diffuse filtering

### 3. Interactive Peak Finding

**Entry point**: `heisenberg run config.yaml -i` (Python) / TUI within IDL

Launches an interactive interface to visually inspect peak detection and tune parameters:
- Adjust `npixmin` (minimum pixels per peak)
- Adjust `nsigma` (significance threshold)
- Modify contour levels (logrange, logspacing)
- Preview peaks before running full analysis

## Command-Line Reference

```
heisenberg run INPUT_FILE [OPTIONS]

Options:
  -n, --no-diffuse-filtering  Skip iterative diffuse filtering
  -i, --interactive           Run interactive peak finding TUI
  -o, --output PATH           Output file path (default: <galaxy>_results.json)
  -d, --output-dir PATH       Output directory for plots and regions
  -p, --plot                  Generate diagnostic plots
  -r, --regions               Export peak positions to DS9 region files
  -v, --verbose               Enable verbose output

heisenberg validate INPUT_FILE    Validate configuration file
heisenberg create-input OUTPUT    Create template configuration file
heisenberg quick-test [OPTIONS]   Run quick test with synthetic data
```

See [CLI Reference](docs/python/cli-reference.md) for complete documentation.

## Configuration

Configuration can be provided in **YAML** or **IDL input_file** format. Key parameter groups:

| Group | Description | Key Parameters |
|-------|-------------|----------------|
| **files** | Input file paths | `datadir`, `starfile`, `gasfile` |
| **basic_map** | Galaxy properties | `distance`, `inclination`, `centrex/y` |
| **aperture** | Aperture settings | `lapmin`, `lapmax`, `naperture` |
| **peak_id** | Peak detection | `npixmin`, `nsigma`, `logrange`, `logspacing` |
| **timeline** | Timescales | `tstariso`, `tgasmini`, `tgasmaxi` |
| **fitting** | Fitting parameters | `nmc`, `ndepth`, `ntry` |
| **iteration** | Diffuse filtering | `iter_criterion`, `iter_nmax` |
| **flags1-4** | Module switches | See full documentation |

See [Configuration Reference](docs/python/configuration.md) for all ~180 parameters with descriptions and typical values.

## Example Configuration Files

- [`Python/examples/example_config.yaml`](Python/examples/example_config.yaml) - YAML format with comments
- [`Python/examples/example_config.inp`](Python/examples/example_config.inp) - IDL input_file format
- [`IDL/input_file`](IDL/input_file) - Original IDL format with full documentation

## Output

The analysis produces several output files:

| Output | Format | Description |
|--------|--------|-------------|
| Results | JSON | All fitted parameters and derived quantities |
| Summary | TXT | Human-readable summary of key results |
| Plots | PNG/PDF | Tuning fork diagram, parameter PDFs, maps |
| Regions | DS9 .reg | Peak positions for visualization |

See [Output Formats](docs/python/output-formats.md) for detailed format specifications.

## Documentation

### Getting Started
- [Installation Guide](docs/installation.md)
- [Quick Start Tutorial](docs/quickstart.md)
- [Scientific Background](docs/scientific-background.md)

### Python Reference
- [CLI Reference](docs/python/cli-reference.md)
- [API Reference](docs/python/api-reference.md)
- [Configuration Reference](docs/python/configuration.md)
- [Interactive Mode Guide](docs/python/interactive-mode.md)

### IDL Reference
- [IDL Usage Guide](docs/idl/usage.md)
- [IDL Configuration](docs/idl/configuration.md)
- [IDL Procedures Reference](docs/idl/procedures.md)

### Tutorials
- [Basic Analysis Tutorial](docs/tutorials/basic-analysis.md)
- [Parameter Tuning Best Practices](docs/tutorials/parameter-tuning.md)
- [Advanced Usage](docs/tutorials/advanced-usage.md)
- [Interpreting Results](docs/tutorials/interpreting-results.md)

### Other
- [Troubleshooting](docs/troubleshooting.md)
- [FAQ](docs/faq.md)
- [Glossary](docs/glossary.md)

## Jupyter Notebook Demo

See [`Python/examples/demo_analysis.ipynb`](Python/examples/demo_analysis.ipynb) for an interactive demonstration including:
- Loading and visualizing FITS data
- Running the full analysis pipeline
- Interpreting and plotting results
- Parameter exploration and sensitivity analysis

## Citation

**Citation is required** when using Heisenberg in any published work, including papers, conference proceedings, theses, and technical reports. Please cite both papers below:

> Proper citation is a condition of use. It acknowledges the significant effort that went into developing both the scientific methodology and this software implementation.

**Method paper:**
```bibtex
@ARTICLE{Kruijssen2014,
    author = {{Kruijssen}, J.~M.~D. and {Longmore}, S.~N.},
    title = "{An uncertainty principle for star formation - I. Why galactic star formation relations break down below a certain spatial scale}",
    journal = {MNRAS},
    year = 2014,
    month = apr,
    volume = 439,
    pages = {3239-3252},
    doi = {10.1093/mnras/stu098},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2014MNRAS.439.3239K}
}
```

**Code paper:**
```bibtex
@ARTICLE{Kruijssen2018,
    author = {{Kruijssen}, J.~M.~D. and {Schruba}, A. and {Hygate}, A.~P.~S. and
              {Hu}, C.-Y. and {Haydon}, D.~T. and {Longmore}, S.~N.},
    title = "{Fast and inefficient star formation due to short-lived molecular clouds and rapid feedback}",
    journal = {MNRAS},
    year = 2018,
    month = sep,
    volume = 479,
    pages = {1866-1952},
    doi = {10.1093/mnras/sty1128},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.1866K}
}
```

### Additional Citations

Please also cite the following papers where applicable:

**If using the Python implementation:**
```bibtex
@MISC{Chevance2026,
    author = {{Chevance}, M. and {Kruijssen}, J.~M.~D.},
    title = "{Heisenberg: Python implementation of the uncertainty principle for star formation}",
    howpublished = {Astrophysics Source Code Library, ascl:XXXX.XXX},
    year = 2026,
    note = {Details TBD}
}
```

**If using reference timescale estimates (tstariso calibration):**
```bibtex
@ARTICLE{Haydon2020,
    author = {{Haydon}, D.~T. and {Kruijssen}, J.~M.~D. and {Chevance}, M. and
              {Hygate}, A.~P.~S. and {Krumholz}, M.~R. and {Schruba}, A. and
              {Longmore}, S.~N.},
    title = "{An uncertainty principle for star formation - III. The characteristic emission time-scales of star formation rate tracers}",
    journal = {MNRAS},
    year = 2020,
    month = oct,
    volume = 498,
    pages = {235-257},
    doi = {10.1093/mnras/staa2430},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2020MNRAS.498..235H}
}
```

**If using iterative diffuse emission filtering:**
```bibtex
@ARTICLE{Hygate2019,
    author = {{Hygate}, A.~P.~S. and {Kruijssen}, J.~M.~D. and {Chevance}, M. and
              {Schruba}, A. and {Haydon}, D.~T. and {Longmore}, S.~N.},
    title = "{An uncertainty principle for star formation - IV. On the nature and filtering of diffuse emission}",
    journal = {MNRAS},
    year = 2019,
    month = sep,
    volume = 488,
    pages = {2800-2824},
    doi = {10.1093/mnras/stz1779},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.2800H}
}
```

## License

Heisenberg is released under the **GNU General Public License v3.0** (GPL-3.0).

This means you are free to use, modify, and distribute this software, provided that:
- Any derivative work is also released under GPL-3.0
- The original copyright and license notices are preserved
- Source code is made available when distributing the software

See [LICENSE](LICENSE) for the full license text.

**Note:** Use of this software in patents requires involvement of the original authors. Please contact us before filing any patent applications that incorporate Heisenberg.

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for:
- Development environment setup
- Code style guidelines
- Testing requirements
- Pull request process

## Authors

**Lead developers:**
- J. M. Diederik Kruijssen
- Mélanie Chevance

**Contributors:**
- Andreas Schruba
- Steven N. Longmore
- Alexander P. S. Hygate
- Daniel T. Haydon
- Jacob L. Ward
- Jaeyeon Kim

## Support

- **Issues**: Please report bugs and feature requests via GitHub Issues
- **Questions**: For usage questions, please use GitHub Discussions

## Acknowledgements

Development of Heisenberg has been supported by ERC-StG-714907 MUSTANG, DFG Emmy Noether grant KR4801/1-1, and DFG Emmy Noether grant CH2137/1-1.

We thank all contributors to the development and testing of this code.

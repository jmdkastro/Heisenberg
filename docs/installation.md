# Installation Guide

This guide covers installation of both the Python (recommended) and IDL (legacy) implementations of Heisenberg.

## Python Installation

### Requirements

- **Python 3.9 or higher**
- Operating System: Linux, macOS, or Windows

### Core Dependencies

The following packages are required and will be installed automatically:

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | >=1.20 | Array operations, FFT |
| scipy | >=1.7 | Optimization, image processing |
| astropy | >=5.0 | FITS I/O, WCS, units |
| matplotlib | >=3.5 | Plotting |
| click | >=8.0 | Command-line interface |
| pyyaml | >=6.0 | YAML configuration parsing |

### Optional Dependencies

| Package | Purpose |
|---------|---------|
| regions | DS9 region file support |
| tqdm | Progress bars |
| rich | Enhanced terminal output for TUI |
| jupyter | Running demo notebooks |

### Installation Methods

#### Method 1: From Source (Recommended for Development)

```bash
# Clone the repository
git clone https://github.com/username/heisenberg.git
cd heisenberg

# Install in development mode
cd Python/
pip install -e .

# Or with all optional dependencies
pip install -e ".[dev]"
```

#### Method 2: Direct pip install

```bash
# If published to PyPI (future)
pip install heisenberg
```

### Verifying the Installation

After installation, verify that everything is working:

```bash
# Check CLI is available
heisenberg --version

# Run quick test with synthetic data
heisenberg quick-test --plot

# Validate example configuration
heisenberg validate ../Python/examples/example_config.yaml
```

In Python:

```python
import heisenberg
print(heisenberg.__version__)

# Test core imports
from heisenberg.core import run_tuningfork, TuningForkConfig
from heisenberg.io import read_fits
print("All imports successful!")
```

### Troubleshooting Python Installation

#### ImportError: No module named 'heisenberg'

Make sure you installed the package:

```bash
pip install -e Python/
```

#### Missing optional dependencies

For full functionality including DS9 region support:

```bash
pip install regions rich tqdm
```

#### Matplotlib backend issues

If plots don't display, try setting the backend:

```python
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'Agg' for non-interactive
```

---

## IDL Installation

### Requirements

- **IDL 8.0 or higher** (Harris Geospatial Solutions)
- **IDL Astronomy User's Library** (required)

### Installing the IDL Astronomy Library

The [IDL Astronomy User's Library](https://idlastro.gsfc.nasa.gov/) provides essential astronomical routines used by Heisenberg.

1. Download from: https://idlastro.gsfc.nasa.gov/

2. Unpack to a directory (e.g., `/usr/local/idl/astrolib/`)

3. Add to your IDL path in your startup file (`~/.idl/idl_startup.pro`):

```idl
!PATH = !PATH + ':' + expand_path('+/usr/local/idl/astrolib/pro/')
```

### Installing Heisenberg IDL

1. Clone or download the Heisenberg repository

2. Add the IDL directory to your IDL path:

```idl
; In your IDL startup file or interactively:
!PATH = !PATH + ':' + expand_path('+/path/to/heisenberg/IDL/')
```

Or set the environment variable before starting IDL:

```bash
export IDL_PATH="${IDL_PATH}:+/path/to/heisenberg/IDL/"
```

### Verifying the IDL Installation

```idl
IDL> ; Test that procedures compile
IDL> .compile tuningfork
IDL> .compile diffuse_iteration
IDL> print, 'Heisenberg IDL installation successful!'
```

### Optional: DS9 for Visualization

For interactive peak visualization, install SAOImage DS9:

- **macOS**: `brew install --cask saoimageds9`
- **Linux**: Available from package managers or https://ds9.si.edu/
- **Windows**: Download from https://ds9.si.edu/

Ensure `ds9` is in your system PATH.

### Troubleshooting IDL Installation

#### Procedure not found errors

Ensure the IDL path is set correctly:

```idl
IDL> print, !PATH
; Should include /path/to/heisenberg/IDL/
```

#### Missing Astronomy Library routines

If you see errors about missing routines like `READFITS` or `SXPAR`:

1. Verify the Astronomy Library is installed
2. Check it's in your IDL_PATH
3. Test: `IDL> .compile readfits`

#### X11 display issues (Linux/macOS)

For GUI features, ensure X11 is properly configured:

```bash
# Test X11
xeyes  # Should open a window

# If using SSH, enable X forwarding
ssh -X user@host
```

---

## Development Installation

For contributing to Heisenberg development:

```bash
# Clone with full history
git clone https://github.com/username/heisenberg.git
cd heisenberg

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows

# Install with development dependencies
cd Python/
pip install -e ".[dev]"

# Install pre-commit hooks (optional)
pre-commit install

# Run tests to verify
pytest tests/
```

### Development Dependencies

The `[dev]` extra installs:

| Package | Purpose |
|---------|---------|
| pytest | Testing framework |
| pytest-cov | Coverage reporting |
| black | Code formatting |
| isort | Import sorting |
| flake8 | Linting |
| mypy | Type checking |

---

## Platform-Specific Notes

### macOS

If using Homebrew Python:

```bash
brew install python@3.11
pip3 install -e Python/
```

### Linux

On Ubuntu/Debian, you may need:

```bash
sudo apt-get install python3-dev python3-pip
```

### Windows

Use Anaconda or Miniconda for easier dependency management:

```bash
conda create -n heisenberg python=3.11
conda activate heisenberg
pip install -e Python/
```

---

## Next Steps

After installation:

1. Read the [Quick Start Guide](quickstart.md) to run your first analysis
2. Explore the [example configurations](../Python/examples/)
3. Try the [demo Jupyter notebook](../Python/examples/demo_analysis.ipynb)

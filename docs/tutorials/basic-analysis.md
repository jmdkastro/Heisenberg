# Basic Analysis Tutorial

Step-by-step guide to running your first Heisenberg analysis.

## Overview

This tutorial walks through a complete analysis from data preparation to interpreting results. We'll use synthetic data for demonstration, but the workflow applies identically to real observations.

## Prerequisites

- Heisenberg installed (see [Installation Guide](../installation.md))
- Two FITS images on the same pixel grid:
  - A star formation tracer (e.g., Hα, FUV, 24μm)
  - A gas tracer (e.g., CO, HI)
- Basic knowledge of your target (distance, centre coordinates)

## Step 1: Verify Installation

First, verify Heisenberg is installed correctly:

```bash
# Check CLI is available
heisenberg --version

# Run quick test with synthetic data
heisenberg quick-test --plot
```

You should see output showing fitted parameters for the synthetic galaxy.

## Step 2: Prepare Your Data

### Data Requirements

Your FITS files must:
1. **Be on the same pixel grid** (same dimensions, same WCS)
2. **Have valid WCS** for pixel-to-physical scale conversion
3. **Trace successive evolutionary phases** (gas → stars)

### Checking Your Data

```python
from astropy.io import fits

# Load and check your files
star_hdu = fits.open('star_map.fits')
gas_hdu = fits.open('gas_map.fits')

print(f"Star map shape: {star_hdu[0].data.shape}")
print(f"Gas map shape: {gas_hdu[0].data.shape}")

# Shapes must match!
assert star_hdu[0].data.shape == gas_hdu[0].data.shape, "Maps must have same dimensions"
```

### Regridding (if needed)

If your maps have different pixel grids, regrid them first using tools like:
- `reproject` (Python/astropy)
- `MONTAGE`
- Your favourite image processing tool

## Step 3: Create Configuration File

### Option A: Copy and Edit Example

```bash
cp Python/examples/example_config.yaml my_analysis.yaml
```

Edit `my_analysis.yaml` with your parameters.

### Option B: Create from Template

```bash
heisenberg create-input my_analysis.yaml
```

### Minimum Required Edits

At minimum, edit these parameters:

```yaml
# File locations
files:
  datadir: /path/to/your/data/    # Absolute path to data directory
  galaxy: MyGalaxy                 # Name for output files
  starfile: my_star_map.fits       # Star formation tracer
  gasfile: my_gas_map.fits         # Gas tracer

# Galaxy properties
basic_map:
  distance: 10000000.0             # Distance in PARSECS (e.g., 10 Mpc = 10e6 pc)
  centrex: 512                     # Galaxy centre X pixel
  centrey: 512                     # Galaxy centre Y pixel
  maxradius: 5000.0                # Maximum radius in pc

# Aperture range
aperture:
  lapmin: 50.0                     # Minimum aperture (pc) - set to ~2× resolution
  lapmax: 3000.0                   # Maximum aperture (pc) - set to ~0.5× FOV

# Stellar tracer timescale
timeline:
  tstariso: 4.0                    # Reference timescale in Myr (4-5 for Hα)
```

## Step 4: Validate Configuration

Before running the full analysis, validate your configuration:

```bash
heisenberg validate my_analysis.yaml
```

Expected output:
```
Configuration is valid!
  Galaxy: MyGalaxy
  Distance: 10000000 pc
  Data directory: /path/to/your/data/
  Star file: my_star_map.fits
  Gas file: my_gas_map.fits
  Apertures: 9 (50 - 3000 pc)
```

If validation fails, check the error message and fix the configuration.

## Step 5: Run Initial Analysis

### First Run: Single-Pass (Quick)

For a quick initial look, run without diffuse filtering:

```bash
heisenberg run my_analysis.yaml --no-diffuse-filtering --verbose --plot
```

This runs faster and gives you initial results to check.

### Check the Output

The run produces:
- `MyGalaxy_results.json` - Full results
- `MyGalaxy_results.txt` - Human-readable summary
- Plots (if `--plot` was used)

Review the summary:
```bash
cat MyGalaxy_results.txt
```

Example output:
```
Heisenberg Analysis Results
===========================
Galaxy: MyGalaxy
Distance: 10.0 Mpc

Fitted Parameters:
  t_gas  = 12.5 (+3.2/-2.8) Myr
  t_over = 2.1 (+0.8/-0.6) Myr
  lambda = 180 (+45/-35) pc

Quality:
  chi2_min = 4.2
  Star peaks: 45
  Gas peaks: 52
```

## Step 6: Verify Peak Identification

The most critical step is verifying that peaks are correctly identified.

### Option A: Check Output Plots

If you ran with `--plot`, examine the peak maps to see detected peaks overlaid on the images.

### Option B: Interactive Mode

For hands-on verification, use interactive mode:

```bash
heisenberg run my_analysis.yaml --interactive --verbose
```

This launches a TUI where you can:
- View peaks on both maps
- Adjust detection parameters
- Accept or reject configurations

### What to Look For

**Good peak identification:**
- All obvious emission peaks are detected
- No spurious noise peaks
- Peaks aren't over-split or over-merged

**If peaks look wrong:**
- Too few peaks → lower `nsigma` (e.g., 5 → 3)
- Too many spurious peaks → raise `nsigma` (e.g., 5 → 7)
- Over-split peaks → increase `logspacing`
- Over-merged peaks → decrease `logspacing`

## Step 7: Run Full Analysis

Once you're satisfied with peak identification, run the full analysis with diffuse filtering:

```bash
heisenberg run my_analysis.yaml --verbose --plot --regions
```

This:
1. Identifies peaks
2. Runs initial tuningfork fit
3. Filters diffuse emission at scale λ
4. Re-fits and iterates until λ converges
5. Calculates derived quantities
6. Generates output files

### Monitor Progress

With `--verbose`, you'll see progress updates:

```
Loading configuration from my_analysis.yaml...
Galaxy: MyGalaxy
Distance: 10000000 pc
Loading stellar map: /path/to/data/my_star_map.fits
Loading gas map: /path/to/data/my_gas_map.fits
Astrometry check passed: star and gas maps match
Image size: (1024, 1024)
Pixel scale: 10.25 pc/pixel

Running iterative diffuse filtering analysis...
Max iterations: 10
Convergence criterion: 5.0%

Iteration 1: lambda = 180 pc
Iteration 2: lambda = 165 pc
Iteration 3: lambda = 168 pc
Iteration 4: lambda = 167 pc
Converged after 4 iterations

============================================================
RESULTS
============================================================
t_gas    = 12.5 (+3.2/-2.8) Myr
t_over   = 2.1 (+0.8/-0.6) Myr
lambda   = 167 (+42/-33) pc
chi2_min = 3.8
Peaks: 45 stellar, 52 gas

Results saved to MyGalaxy_results.json
Summary saved to MyGalaxy_results.txt
Plots saved to . (3 files)
Region files saved to . (2 files)
```

## Step 8: Interpret Results

### Primary Fitted Parameters

| Parameter | Symbol | Meaning |
|-----------|--------|---------|
| `tgas` | t_gas | Gas cloud lifetime from assembly to dispersal |
| `tover` | t_over | Overlap time (embedded phase duration) |
| `lambda_` | λ | Region separation length |
| `chi2_min` | χ² | Goodness of fit (lower is better) |

### Quality Assessment

**Good fit indicators:**
- χ² close to 1 (typically 1-10 is acceptable)
- Smooth tuning fork diagram
- Unimodal parameter PDFs
- λ > 1.5 × lapmin

**Warning signs:**
- Very high χ² (>50): poor model fit
- Bimodal PDFs: ambiguous solution
- λ ≈ lapmin: resolution limited
- Very large error bars: insufficient statistics

### Derived Quantities

If conversion factors are provided, additional quantities are derived:

| Quantity | Description |
|----------|-------------|
| `esf` | Star formation efficiency per free-fall time |
| `vfb` | Feedback velocity |
| `mass_loading` | Mass loading factor |

## Step 9: Export and Visualize

### DS9 Region Files

If you ran with `--regions`, peak positions are exported:
- `MyGalaxy_star_peaks.reg`
- `MyGalaxy_gas_peaks.reg`

View in DS9:
```bash
ds9 my_star_map.fits -region MyGalaxy_star_peaks.reg
```

### Reading Results in Python

```python
import json

with open('MyGalaxy_results.json') as f:
    results = json.load(f)

# Access fitted parameters
print(f"t_gas = {results['fit']['tgas']:.1f} Myr")
print(f"t_over = {results['fit']['tover']:.1f} Myr")
print(f"lambda = {results['fit']['lambda_']:.0f} pc")

# Access errors
print(f"t_gas error: +{results['fit']['tgas_errmax']:.1f}/-{results['fit']['tgas_errmin']:.1f}")

# Access peak counts
print(f"Star peaks: {len(results['star_peaks'])}")
print(f"Gas peaks: {len(results['gas_peaks'])}")
```

## Step 10: Iterate if Needed

Based on results, you may need to adjust parameters and re-run:

| Issue | Solution |
|-------|----------|
| Poor peak identification | Adjust `npixmin`, `nsigma`, `logspacing` |
| λ too close to lapmin | Decrease `lapmin` or increase `peak_res` |
| Non-converging iterations | Increase `iter_criterion` or check data |
| High χ² | Check data quality, aperture range |

See [Parameter Tuning](parameter-tuning.md) for detailed guidance.

## Complete Example Script

```bash
#!/bin/bash
# Complete Heisenberg analysis workflow

GALAXY="NGC1234"
CONFIG="${GALAXY}_config.yaml"
OUTDIR="results/${GALAXY}"

# Create output directory
mkdir -p "$OUTDIR"

# 1. Validate configuration
echo "Validating configuration..."
heisenberg validate "$CONFIG"

# 2. Quick initial run (no diffuse filtering)
echo "Running initial analysis..."
heisenberg run "$CONFIG" \
    --no-diffuse-filtering \
    --output "${OUTDIR}/${GALAXY}_initial.json" \
    --verbose

# 3. Full analysis with diffuse filtering
echo "Running full analysis..."
heisenberg run "$CONFIG" \
    --output "${OUTDIR}/${GALAXY}_results.json" \
    --output-dir "$OUTDIR" \
    --plot \
    --regions \
    --verbose

echo "Analysis complete! Results in ${OUTDIR}/"
```

## Next Steps

- [Parameter Tuning](parameter-tuning.md) - Best practices for parameter selection
- [Interpreting Results](interpreting-results.md) - Understanding output parameters
- [Advanced Usage](advanced-usage.md) - Diffuse filtering, masking, batch processing
- [Demo Notebook](../../Python/examples/demo_analysis.ipynb) - Interactive exploration

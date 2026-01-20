# Quick Start Guide

Get up and running with Heisenberg in 5 minutes.

## Prerequisites

- Heisenberg installed (see [Installation Guide](installation.md))
- Two FITS images: a star formation tracer and a gas tracer for the same galaxy
- Basic knowledge of your data (distance, pixel scale, tracer timescales)

## 1. Quick Test with Synthetic Data

First, verify your installation by running a quick test with synthetic data:

```bash
heisenberg quick-test --plot
```

This generates synthetic star and gas maps, runs the tuningfork analysis, and displays the results. You should see output like:

```
Generating synthetic galaxy with 15 peaks...
Running tuningfork analysis...

==================================================
RESULTS
==================================================
t_gas    = 8.52 (+2.15/-1.83) Myr
t_over   = 1.23 (+0.45/-0.38) Myr
lambda   = 45.2 (+8.1/-7.3) pc
chi2_min = 3.45
Peaks: 12 stellar, 14 gas
```

## 2. Prepare Your Configuration File

Copy an example configuration file and modify it for your data:

```bash
cp Python/examples/example_config.yaml my_config.yaml
```

Edit the key parameters:

```yaml
# File paths
files:
  datadir: /path/to/your/data/
  galaxy: NGC1234
  starfile: ngc1234_halpha.fits   # Star formation tracer
  gasfile: ngc1234_co.fits        # Gas tracer

# Galaxy properties
basic_map:
  distance: 10000000.0   # Distance in parsecs (10 Mpc)
  inclination: 30.0      # Inclination in degrees
  centrex: 512           # Galaxy center X pixel
  centrey: 512           # Galaxy center Y pixel

# Aperture settings
aperture:
  lapmin: 50.0           # Minimum aperture diameter (pc)
  lapmax: 3200.0         # Maximum aperture diameter (pc)
  naperture: 9           # Number of aperture sizes

# Timeline - IMPORTANT: set for your stellar tracer
timeline:
  tstariso: 4.0          # Stellar tracer timescale (Myr)
                         # H-alpha: ~4-5 Myr
                         # FUV: ~100 Myr
                         # 24um: ~10-20 Myr
```

## 3. Validate Your Configuration

Before running the analysis, validate your configuration:

```bash
heisenberg validate my_config.yaml
```

You should see:

```
Configuration is valid!
  Galaxy: NGC1234
  Distance: 10000000 pc
  Data directory: /path/to/your/data/
  Star file: ngc1234_halpha.fits
  Gas file: ngc1234_co.fits
  Apertures: 9 (50 - 3200 pc)
```

## 4. Run the Analysis

### Option A: Full Analysis with Diffuse Filtering (Recommended)

```bash
heisenberg run my_config.yaml --verbose --plot
```

This runs iterative diffuse filtering, which:
1. Identifies peaks in both maps
2. Fits the KL14 model to get initial lambda
3. Filters diffuse emission at scale lambda
4. Refits and repeats until lambda converges

### Option B: Single-Pass Analysis (Quick)

For a quick initial look without diffuse filtering:

```bash
heisenberg run my_config.yaml --no-diffuse-filtering --verbose --plot
```

### Option C: Interactive Peak Finding

If you want to tune peak detection parameters:

```bash
heisenberg run my_config.yaml --interactive
```

This launches an interactive interface where you can:
- View detected peaks overlaid on maps
- Adjust detection thresholds
- Preview results before committing

## 5. Understand the Output

The analysis produces several files:

| File | Contents |
|------|----------|
| `NGC1234_results.json` | All fitted and derived parameters |
| `NGC1234_results.txt` | Human-readable summary |
| `NGC1234_tuningfork.png` | Tuning fork diagram (if --plot) |

### Key Output Parameters

| Parameter | Symbol | Meaning |
|-----------|--------|---------|
| `tgas` | t_gas | Gas cloud lifetime (Myr) |
| `tover` | t_over | Overlap/feedback time (Myr) |
| `lambda_` | λ | Region separation length (pc) |
| `chi2_min` | χ² | Goodness of fit |

### Reading Results in Python

```python
import json

with open('NGC1234_results.json') as f:
    results = json.load(f)

print(f"Gas lifetime: {results['fit']['tgas']:.1f} Myr")
print(f"Overlap time: {results['fit']['tover']:.1f} Myr")
print(f"Lambda: {results['fit']['lambda_']:.0f} pc")
```

## 6. Running in IDL

If you prefer IDL:

```idl
; With diffuse filtering (recommended)
IDL> inputfile = '/path/to/my_config.inp'
IDL> .run heisenberg

; Without diffuse filtering
IDL> inputfile = '/path/to/my_config.inp'
IDL> .run heisenberg_nodf
```

## Common Issues

### "Too few peaks found"

- Lower `nsigma` (significance threshold) from 5 to 3-4
- Lower `npixmin` (minimum pixels per peak) from 20 to 10-15
- Check that masks aren't excluding too much area

### "Non-converging iterations"

- Check data quality (high noise, artifacts)
- Try different `initial_guess` for lambda
- Increase `iter_criterion` from 0.05 to 0.10

### "Results seem unreasonable"

- Verify `tstariso` is correct for your stellar tracer
- Check `distance` is in parsecs (not Mpc!)
- Ensure `lapmin` > 2× your resolution

## Next Steps

- Read the [Basic Analysis Tutorial](tutorials/basic-analysis.md) for a detailed walkthrough
- Learn about [Parameter Tuning](tutorials/parameter-tuning.md) best practices
- Explore the [Configuration Reference](python/configuration.md) for all parameters
- Try the [Demo Notebook](../Python/examples/demo_analysis.ipynb) for interactive exploration

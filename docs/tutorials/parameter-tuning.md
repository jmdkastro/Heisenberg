# Parameter Tuning Best Practices

This guide provides recommended parameter values and strategies for achieving convergent, reliable fits with Heisenberg. Based on extensive experience running the code on diverse datasets.

## Overview

Most parameters should be left at their default values. This guide highlights the parameters that typically need adjustment for each new analysis and provides guidance on selecting appropriate values.

## Parameters to Adapt for Each Galaxy

### Galaxy Properties (INPUT PARAMETERS 1)

| Parameter | Description | How to Set |
|-----------|-------------|------------|
| `distance` | Distance to galaxy in pc | From literature or your measurements |
| `inclination` | Inclination angle in degrees | From literature or morphological fitting |
| `posangle` | Position angle in degrees | From literature or morphological fitting |
| `centrex`, `centrey` | Galaxy centre pixel coordinates | From image inspection or header |
| `maxradius` | Maximum analysis radius in pc | Approximate total field of view (not much larger) |

**Example:**
```yaml
basic_map:
  distance: 840000.0      # M31 distance in pc
  inclination: 77.0       # M31 inclination
  posangle: 38.0          # M31 position angle
  centrex: 701
  centrey: 701
  minradius: 0.0
  maxradius: 10000.0      # Approximate field of view
```

### Aperture Settings (INPUT PARAMETERS 2)

| Parameter | Recommendation |
|-----------|----------------|
| `lapmin` | Set to the **worst spatial resolution** of your two maps |
| `lapmax` | Set to **approximate total field of view** (not much more) |
| `peak_res` | Set to **0** (use smallest aperture for peak identification) |
| `max_res` | Set to **naperture - 1** (use largest aperture for galactic averages) |

**Critical rule:** λ should be at least **1.5× the smallest aperture size** (`lapmin`). If fitted λ is close to `lapmin`, either decrease `lapmin` or increase `peak_res`.

### Peak Identification (INPUT PARAMETERS 3)

These parameters require case-by-case tuning. **Always verify peak identification** by either:
1. Activating the interactive mode (`peak_find_tui = 1`)
2. Checking the output figures: `figures/***_map_gas.eps` and `figures/***_map_star.eps`

| Parameter | Default | Typical Range | Description |
|-----------|---------|---------------|-------------|
| `npixmin` | 20 | 10-50 | Minimum pixels per peak |
| `nsigma` | 5 | 3-10 | Significance threshold |
| `logrange_s` | 2.0 | 1.5-3.0 | Log range for stellar contours |
| `logspacing_s` | 0.5 | 0.1-0.5 | Log spacing for stellar contours |
| `logrange_g` | 2.0 | 1.5-3.0 | Log range for gas contours |
| `logspacing_g` | 0.5 | 0.1-0.5 | Log spacing for gas contours |

**Tuning strategy:**
1. Start with defaults
2. Run analysis and check peak maps
3. If too few peaks: lower `nsigma` or `npixmin`
4. If too many spurious peaks: raise `nsigma` or `npixmin`
5. If peaks are over-split: increase `logspacing`
6. If peaks are over-merged: decrease `logspacing`

### Timeline Parameters (INPUT PARAMETERS 4)

| Parameter | Recommendation |
|-----------|----------------|
| `tstariso` | Set according to your stellar tracer (see table below) |
| `tgasmaxi` | Can be lowered to ~200 Myr (typically sufficient) |

**Reference timescales by tracer** (from Haydon et al. 2020):

| Stellar Tracer | tstariso (Myr) | Notes |
|----------------|----------------|-------|
| Hα | 4-5 | Most common choice |
| FUV | ~100 | Long timescale, constrain tover carefully |
| NUV | 100-200 | Similar to FUV |
| 24 μm | 10-20 | Embedded SF tracer |
| Hα + 24μm | 5-10 | Composite tracer |

### Conversion Factors (INPUT PARAMETERS 7)

| Parameter | Recommendation |
|-----------|----------------|
| `convstar` | Adapt to your galaxy, or leave at 0 if `map_units = 0` |
| `convgas` | Adapt to your galaxy, or leave at 0 if `map_units = 0` |

If you don't have calibrated conversion factors, set `map_units = 0` in FLAGS 4 to work with relative quantities.

---

## Noise Threshold Parameters (INPUT PARAMETERS 9)

These are important for iterative diffuse filtering. **Ignore if running `heisenberg_nodf`.**

| Parameter | Default | Recommendation |
|-----------|---------|----------------|
| `use_noisecut` | 1 | Keep at 1 |
| `noisethresh_s` | 20.0 | Start with 0, then set to ~3× background σ |
| `noisethresh_g` | 20.0 | Start with 0, then set to ~3× background σ |

**Strategy:**
1. Start by setting `noisethresh_s` and `noisethresh_g` to **0**
   - This sets negative pixels after filtering to 0
2. If significant noise/background remains, measure the standard deviation of the background (emission-free regions)
3. Set thresholds to approximately **3× the background standard deviation**

---

## Flags to Change from Defaults

Most flags should remain at their default values. Change only when necessary:

### FLAGS 1 (Module Switches)

| Flag | Default | When to Change |
|------|---------|----------------|
| `mask_images` | 1 | Set to **0** if masks are not provided |
| `cleanup` | 2 | Set to **0** to keep all temporary files (debugging) |

### FLAGS 3 (Masking)

| Flag | Default | When to Change |
|------|---------|----------------|
| `mstar_ext`, `mstar_int` | 1 | Set to **0** if masks are not provided |
| `mgas_ext`, `mgas_int` | 1 | Set to **0** if masks are not provided |
| `cut_radius` | 1 | Set to **0** to keep the entire field of view |

### FLAGS 4 (Analysis Options)

| Flag | Default | When to Change |
|------|---------|----------------|
| `set_centre` | 1 | Set to **0** to use central pixel as galaxy centre |
| `peak_find_tui` | 0 | Set to **1** for interactive peak selection (requires DS9) |
| `map_units` | 1 | Set to **0** if conversion factors are not provided |

---

## Output Verification Checklist

After running Heisenberg, check the following outputs to verify the fit quality:

### 1. Peak Maps
**Files:** `figures/***_map_gas.eps`, `figures/***_map_star.eps`

**Check that:**
- All obvious peaks are identified
- No spurious noise peaks are included
- Peak boundaries are reasonable

**If issues:** Adjust `npixmin`, `nsigma`, `logrange`, `logspacing`

### 2. Tuning Fork Diagram
**File:** `figures/***_fit.eps`

**Check that:**
- The tuning fork is smooth without additional wiggles
- Data points follow the expected tuning fork shape
- Error bars are reasonable
- The fit passes through the data

**If wiggles present:** Adjust diffuse emission filtering parameters

### 3. Parameter PDFs
**Files:** `figures/***_distr_tgas.eps`, `***_distr_tover.eps`, `***_distr_lambda.eps`

**Check that:**
- 1D marginalised PDFs are smooth
- PDFs are unimodal (single peak)
- Parameters are not hitting boundaries

**If bimodal or hitting boundaries:** Expand fitting range or check data quality

### 4. Lambda vs Aperture Scale
**Critical check:** λ should be at least **1.5× lapmin**

**If λ ≈ lapmin:**
- Decrease `lapmin` (if resolution allows)
- Or increase `peak_res` to start fitting from a larger aperture

---

## Iteration Parameters (Diffuse Filtering)

For iterative diffuse filtering (`heisenberg.pro`), these parameters control convergence:

| Parameter | Default | Recommendation |
|-----------|---------|----------------|
| `use_guess` | 0 | Usually keep at 0 |
| `initial_guess` | 200 | Rough estimate of λ if using `use_guess = 1` |
| `iter_criterion` | 0.05 | 5% convergence threshold; increase if oscillating |
| `iter_nmax` | 10 | Maximum iterations; usually converges in 3-7 |

**Typical convergence:** Most analyses converge within 3-7 iterations.

---

## Common Issues and Solutions

### Issue: Too few peaks detected

**Symptoms:** Very few peaks in output maps, poor statistics

**Solutions:**
1. Lower `nsigma` from 5 to 3-4
2. Lower `npixmin` from 20 to 10-15
3. Decrease `logspacing` to create more contour levels
4. Check mask coverage isn't too aggressive

### Issue: Too many spurious peaks

**Symptoms:** Noise peaks identified, scattered peaks in low-emission regions

**Solutions:**
1. Increase `nsigma` from 5 to 7-10
2. Increase `npixmin` from 20 to 30-50
3. Increase `logspacing` to reduce contour levels
4. Add masks to exclude noisy regions

### Issue: Tuning fork has wiggles

**Symptoms:** Non-smooth tuning fork with oscillations

**Solutions:**
1. Adjust diffuse emission filtering
2. Check noise thresholds (`noisethresh_s`, `noisethresh_g`)
3. Verify peak identification is clean

### Issue: λ at or near lapmin

**Symptoms:** Fitted λ very close to minimum aperture size

**Solutions:**
1. Decrease `lapmin` if resolution allows
2. Increase `peak_res` to start fitting from larger aperture
3. Check if data truly has small-scale structure

### Issue: Non-converging iterations

**Symptoms:** λ oscillates without settling

**Solutions:**
1. Increase `iter_criterion` from 0.05 to 0.10
2. Use `use_guess = 1` with better `initial_guess`
3. Check for data artifacts or edge effects

### Issue: Parameters hitting boundaries

**Symptoms:** Fitted tgas or tover at tgasmini/tgasmaxi

**Solutions:**
1. Expand fitting range (tgasmini, tgasmaxi)
2. Check `tstariso` is appropriate for your tracer
3. Verify data quality and calibration

---

## Example Configurations

### Nearby Spiral Galaxy (Hα + CO)

```yaml
files:
  galaxy: NGC628
  starfile: ngc628_halpha.fits
  gasfile: ngc628_co21.fits

basic_map:
  distance: 9800000.0     # 9.8 Mpc
  inclination: 7.0
  maxradius: 12000.0

aperture:
  lapmin: 50.0            # ~2× resolution
  lapmax: 6400.0
  naperture: 9
  peak_res: 0
  max_res: 8

peak_id:
  npixmin: 20
  nsigma: 5.0

timeline:
  tstariso: 4.0           # Hα tracer
```

### Distant Galaxy (lower resolution)

```yaml
basic_map:
  distance: 50000000.0    # 50 Mpc

aperture:
  lapmin: 200.0           # Limited by resolution
  lapmax: 10000.0
  naperture: 7

peak_id:
  npixmin: 10             # Fewer pixels per beam
  nsigma: 4.0             # Lower threshold for S/N
```

### High-Resolution Data (HST + ALMA)

```yaml
aperture:
  lapmin: 20.0            # Excellent resolution
  lapmax: 2000.0
  naperture: 12           # More apertures possible

peak_id:
  npixmin: 30             # More pixels available
  nsigma: 6.0             # Higher threshold possible
```

---

## Summary: Quick Reference

### Parameters to Always Check

1. **Galaxy properties:** `distance`, `inclination`, `centrex/y`, `maxradius`
2. **Apertures:** `lapmin` (set to resolution), `lapmax` (set to FOV)
3. **Tracer timescale:** `tstariso` (from literature)
4. **Peak identification:** Verify with maps or interactive mode

### Parameters to Leave at Defaults

- Most flags (unless masks not provided)
- Fitting parameters (`nmc`, `ndepth`, `ntry`)
- Fourier filter parameters (for `heisenberg_nodf`)
- Most iteration parameters

### Verification Steps

1. Check peak maps for completeness
2. Check tuning fork for smoothness
3. Check PDFs for unimodality
4. Verify λ > 1.5× lapmin

---

## References

- Kruijssen & Longmore (2014), MNRAS 439, 3239
- Kruijssen et al. (2018), MNRAS 479, 1866
- Haydon et al. (2020), MNRAS 498, 235 (tracer timescales)

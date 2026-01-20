# Frequently Asked Questions

## General Questions

### What is Heisenberg?

Heisenberg is a tool that implements the "uncertainty principle for star formation" methodology from Kruijssen & Longmore (2014). It measures the physics of star formation and feedback by analysing the spatial decorrelation between gas and young stellar tracers.

### What can Heisenberg measure?

- **t_gas**: Gas cloud lifetime (from assembly to dispersal)
- **t_over**: Overlap/feedback timescale (embedded phase duration)
- **λ (lambda)**: Characteristic separation between independent star-forming regions
- **Derived quantities**: Star formation efficiency, feedback velocity, mass loading factor

### What data do I need?

Two images on the same pixel grid:
1. A **star formation tracer** (e.g., Hα, FUV, 24μm)
2. A **gas tracer** (e.g., CO, HI)

The images must trace successive phases of the evolutionary sequence (gas → stars).

### Which implementation should I use: Python or IDL?

**Python is recommended** for new users. It has:
- Modern, well-tested codebase
- Better error messages
- Easier installation
- Active development

The **IDL version** is provided for backward compatibility and users with existing IDL workflows.

---

## Data Requirements

### Do my images need to be on the same pixel grid?

Yes. Both images must have:
- Same dimensions (number of pixels)
- Same WCS (same pixel scale and astrometry)

If they differ, regrid one to match the other before running Heisenberg.

### What pixel scale should I use?

There's no strict requirement, but:
- **Minimum**: Apertures should contain multiple pixels
- **Maximum**: Not so large that individual regions are unresolved
- **Typical**: 1-50 pc/pixel for nearby galaxies

### Do I need to mask my images?

Masking is optional but recommended for:
- Excluding foreground stars or background galaxies
- Removing edge effects
- Focusing on specific regions of interest

### What format should my images be in?

FITS format with valid WCS headers.

---

## Parameter Selection

### How do I choose lapmin and lapmax?

- **lapmin**: Set to approximately 2× your spatial resolution (beam FWHM)
- **lapmax**: Set to approximately half your field of view

### What value should I use for tstariso?

This depends on your stellar tracer:

| Tracer | tstariso (Myr) |
|--------|----------------|
| Hα | 4-5 |
| FUV | ~100 |
| 24 μm | 10-20 |

See Haydon et al. (2020) for detailed tracer timescales.

### How do I know if peak identification is correct?

1. Use `--interactive` mode to visually inspect peaks
2. Check the output peak maps
3. Ensure all obvious emission peaks are detected
4. Verify no spurious noise peaks are included

### What if I don't have conversion factors?

Set `map_units = 0` in the configuration. Heisenberg will work with relative quantities (timescale ratios) instead of absolute values.

---

## Running the Analysis

### Should I use diffuse filtering or not?

**With diffuse filtering** (`heisenberg run`):
- Recommended for most scientific analyses
- Properly accounts for extended emission
- Iterates until λ converges

**Without diffuse filtering** (`heisenberg run --no-diffuse-filtering`):
- Faster, good for initial exploration
- Useful for pre-filtered data
- Single-pass analysis

### How long does the analysis take?

Typical runtimes:
- Quick test (synthetic data): ~30 seconds
- Single-pass (no filtering): 1-5 minutes
- Full iterative analysis: 5-30 minutes

Factors affecting runtime:
- Image size
- Number of Monte Carlo samples (`nmc`)
- Number of iterations
- Grid resolution (`ntry`)

### Why is my analysis not converging?

Common causes:
1. Poor initial guess for λ
2. Noisy data creating unstable peaks
3. Inappropriate aperture range
4. Edge effects or artifacts

Solutions:
- Increase `iter_criterion` (e.g., 0.05 → 0.10)
- Provide `initial_guess` with `use_guess = true`
- Check data quality and masks

---

## Interpreting Results

### What is a good chi-squared value?

- **χ² ≈ 1-10**: Good fit
- **χ² > 50**: Poor fit, check data/parameters
- **χ² < 1**: Possible overfitting or error underestimation

### What if λ is close to lapmin?

This suggests the fitted scale is near your resolution limit:
1. Decrease `lapmin` if resolution allows
2. Increase `peak_res` to exclude the smallest apertures from fitting
3. Consider if your data truly has small-scale structure

### What if parameters hit the fitting boundaries?

Expand the fitting range:
```yaml
timeline:
  tgasmini: 0.01   # Lower minimum
  tgasmaxi: 500.0  # Higher maximum
```

### How do I assess the reliability of my results?

Check:
1. **Parameter PDFs**: Should be smooth and unimodal
2. **Tuning fork**: Should be smooth without wiggles
3. **Error bars**: Reasonable (not too large or small)
4. **Peak maps**: All obvious peaks detected
5. **χ²**: Close to 1

---

## Comparison with Other Methods

### How does Heisenberg compare to cloud-by-cloud matching?

Heisenberg:
- Uses statistical approach (many regions)
- Doesn't require individual cloud identification
- Works on unresolved or marginally resolved data
- Provides population-averaged timescales

Cloud-by-cloud matching:
- Requires resolved individual clouds
- Provides cloud-specific measurements
- More sensitive to completeness issues

### Can I use Heisenberg on simulation data?

Yes! Heisenberg works on any two maps tracing successive phases. For simulations:
- Create synthetic observables from simulation output
- Apply appropriate noise and resolution effects
- Compare with known input timescales

---

## Technical Questions

### Why do I get different results with different random seeds?

Monte Carlo sampling introduces statistical variation. For reproducible results:
- Set a fixed seed in configuration
- Increase `nmc` for better statistics
- Report uncertainties from multiple runs

### Can I run Heisenberg in parallel?

The current implementation is single-threaded. For batch processing:
- Run multiple analyses sequentially
- Or use shell-level parallelization (different processes for different galaxies)

### How do I cite Heisenberg?

Please cite:

1. **Method paper**: Kruijssen & Longmore (2014), MNRAS 439, 3239
2. **Code paper**: Kruijssen et al. (2018), MNRAS 479, 1866

---

## Getting Help

### Where can I report bugs?

Open an issue on GitHub with:
- Heisenberg version
- Python/IDL version
- Operating system
- Complete error message
- Minimal example to reproduce

### Where can I ask questions?

- Check this FAQ first
- Search existing GitHub issues
- Open a new discussion on GitHub

### Is there a mailing list or forum?

Currently, GitHub Discussions is the primary forum for questions and discussions.

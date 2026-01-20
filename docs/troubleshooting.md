# Troubleshooting

Common issues and solutions when running Heisenberg.

## Installation Issues

### Python: "ModuleNotFoundError: No module named 'heisenberg'"

**Cause:** Package not installed or not in Python path.

**Solution:**
```bash
cd Python/
pip install -e .

# Verify installation
python -c "import heisenberg; print(heisenberg.__version__)"
```

### Python: "command not found: heisenberg"

**Cause:** CLI script not in PATH.

**Solution:**
```bash
# Find where pip installs scripts
python -m site --user-base
# Add the bin directory to PATH, e.g.:
export PATH="$HOME/.local/bin:$PATH"

# Or run via Python
python -m heisenberg.cli --help
```

### IDL: "% Procedure not found"

**Cause:** IDL path not configured correctly.

**Solution:**
```idl
; Add Heisenberg to path
!PATH = !PATH + ':' + expand_path('+/path/to/heisenberg/IDL/')

; Verify
.compile tuningfork
```

### IDL: Missing Astronomy Library routines

**Cause:** IDL Astronomy User's Library not installed.

**Solution:**
1. Download from https://idlastro.gsfc.nasa.gov/
2. Add to IDL_PATH:
```idl
!PATH = !PATH + ':' + expand_path('+/path/to/astrolib/pro/')
```

---

## Data Issues

### "Astrometry mismatch" / "WCS mismatch"

**Cause:** Star and gas maps have different pixel grids or WCS.

**Solution:**
1. Verify both images have the same dimensions:
```python
from astropy.io import fits
star = fits.open('star.fits')[0].data
gas = fits.open('gas.fits')[0].data
print(star.shape, gas.shape)  # Must match!
```

2. If different, regrid one to match the other:
```python
from reproject import reproject_interp
# Regrid gas to match star
gas_reproj, _ = reproject_interp(gas_hdu, star_header)
```

### "No WCS found"

**Cause:** FITS file lacks WCS information in header.

**Solution:**
1. Add WCS to header if you know the astrometry
2. Or provide pixel scale explicitly in configuration

### "FITS file not found"

**Cause:** Incorrect path or filename.

**Solution:**
1. Use absolute paths in `datadir`
2. Verify filenames match exactly (case-sensitive on Linux/macOS)
3. Check file permissions

---

## Configuration Issues

### "Configuration error: Missing required field"

**Cause:** Required parameter not specified.

**Solution:** Check which field is missing and add it:
```yaml
# These are always required:
files:
  datadir: /path/to/data/
  galaxy: MyGalaxy
  starfile: star.fits
  gasfile: gas.fits

basic_map:
  distance: 10000000.0
  centrex: 512
  centrey: 512
  maxradius: 5000.0
```

### "peak_res must be < naperture"

**Cause:** `peak_res` index exceeds number of apertures.

**Solution:**
```yaml
aperture:
  naperture: 9
  peak_res: 0      # Must be 0 to naperture-1
  max_res: 8       # Must be >= peak_res
```

### "lapmin must be less than lapmax"

**Cause:** Aperture range is inverted.

**Solution:**
```yaml
aperture:
  lapmin: 50.0     # Smaller value
  lapmax: 3000.0   # Larger value
```

---

## Analysis Issues

### Too few peaks detected

**Symptoms:** Very few peaks found, poor statistics.

**Solutions:**
1. Lower significance threshold:
```yaml
peak_id:
  nsigma: 3.0      # Down from 5.0
```

2. Lower minimum pixel count:
```yaml
peak_id:
  npixmin: 10      # Down from 20
```

3. Check masks aren't too aggressive

4. Verify data quality and signal-to-noise

### Too many spurious peaks

**Symptoms:** Noise peaks identified, scattered weak peaks.

**Solutions:**
1. Raise significance threshold:
```yaml
peak_id:
  nsigma: 7.0      # Up from 5.0
```

2. Raise minimum pixel count:
```yaml
peak_id:
  npixmin: 30      # Up from 20
```

3. Add masks to exclude noisy regions

4. Increase contour spacing to reduce over-splitting:
```yaml
peak_id:
  logspacing_s: 0.7  # Up from 0.5
  logspacing_g: 0.7
```

### λ at or near lapmin

**Symptoms:** Fitted λ very close to minimum aperture size.

**Solutions:**
1. Decrease `lapmin` if resolution allows
2. Increase `peak_res` to start fitting from larger aperture:
```yaml
aperture:
  peak_res: 1      # Up from 0
```
3. Verify your resolution estimate is correct

### Non-converging iterations

**Symptoms:** λ oscillates without settling, maximum iterations reached.

**Solutions:**
1. Relax convergence criterion:
```yaml
iteration:
  iter_criterion: 0.10  # Up from 0.05
```

2. Increase maximum iterations:
```yaml
iteration:
  iter_nmax: 20    # Up from 10
```

3. Provide better initial guess:
```yaml
iteration:
  use_guess: true
  initial_guess: 150.0  # Your estimate
```

4. Check data for artifacts (edges, bright point sources)

### Very high chi-squared (χ² > 50)

**Symptoms:** Poor model fit indicated by high χ².

**Solutions:**
1. Check aperture range covers relevant scales
2. Verify peak identification is correct
3. Check for data quality issues
4. Consider if the assumption of independent regions holds

### Parameters hitting boundaries

**Symptoms:** Fitted tgas or tover at tgasmini/tgasmaxi limits.

**Solutions:**
1. Expand fitting range:
```yaml
timeline:
  tgasmini: 0.01   # Lower minimum
  tgasmaxi: 500.0  # Higher maximum
```

2. Verify `tstariso` is appropriate for your tracer

3. Check data calibration

---

## Runtime Issues

### Out of memory

**Cause:** Large images or too many Monte Carlo samples.

**Solutions:**
1. Reduce image size by trimming or rebinning
2. Reduce Monte Carlo samples for testing:
```yaml
fitting:
  nmc: 100         # Down from 1000
```
3. Process smaller regions using `minradius`/`maxradius`

### Analysis takes too long

**Cause:** Large parameter space or many Monte Carlo samples.

**Solutions:**
1. Reduce Monte Carlo samples:
```yaml
fitting:
  nmc: 500         # Moderate value for testing
```

2. Reduce grid resolution:
```yaml
fitting:
  ntry: 51         # Down from 101
```

3. Reduce refinement depth:
```yaml
fitting:
  ndepth: 3        # Down from 4
```

4. Run single-pass first (no diffuse filtering):
```bash
heisenberg run config.yaml --no-diffuse-filtering
```

### X11 display errors (remote/headless)

**Cause:** No display available for plotting.

**Solutions:**
1. Use non-interactive backend:
```python
import matplotlib
matplotlib.use('Agg')
```

2. Disable X11 in IDL configuration:
```
0                    ; use_X11 in FLAGS 4
```

3. Use X11 forwarding:
```bash
ssh -X user@remote
```

---

## Output Issues

### Empty or missing output files

**Cause:** Analysis failed or output disabled.

**Solutions:**
1. Check for error messages in terminal output
2. Verify `write_output` flag is enabled:
```yaml
flags1:
  write_output: true
```
3. Check output directory permissions

### Plots not generated

**Cause:** Matplotlib not available or plotting disabled.

**Solutions:**
1. Install matplotlib:
```bash
pip install matplotlib
```

2. Use `--plot` flag:
```bash
heisenberg run config.yaml --plot
```

3. In IDL, ensure `generate_plot = 1` in FLAGS 1

### Region files not created

**Cause:** Region export not requested.

**Solution:** Use `--regions` flag:
```bash
heisenberg run config.yaml --regions
```

---

## Getting Help

If you can't resolve an issue:

1. **Check the documentation:**
   - [Configuration Reference](python/configuration.md)
   - [Parameter Tuning](tutorials/parameter-tuning.md)
   - [FAQ](faq.md)

2. **Search existing issues:**
   - GitHub Issues

3. **Open a new issue** with:
   - Heisenberg version
   - Python/IDL version
   - Operating system
   - Complete error message
   - Minimal configuration to reproduce

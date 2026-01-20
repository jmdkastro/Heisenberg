# Configuration Reference

Complete reference for all Heisenberg configuration parameters (~180 parameters organized into groups).

## Overview

Configuration can be provided in **YAML format** (recommended) or **IDL input_file format** (for backward compatibility). Parameters are organized into logical groups:

| Group | Description |
|-------|-------------|
| [files](#files) | Input file paths |
| [peak_files](#peak_files) | Peak identification file paths |
| [mask_files](#mask_files) | DS9 mask region files |
| [unfilt_files](#unfilt_files) | Unfiltered images for diffuse calculation |
| [flags1](#flags1) | Module switches |
| [flags2](#flags2) | Ancillary file usage |
| [flags3](#flags3) | Masking options |
| [flags4](#flags4) | Analysis options |
| [basic_map](#basic_map) | Galaxy properties |
| [aperture](#aperture) | Aperture settings |
| [peak_id](#peak_id) | Peak detection parameters |
| [timeline](#timeline) | Timescale parameters |
| [fitting](#fitting) | Fitting parameters |
| [fourier](#fourier) | Fourier filtering parameters |
| [conversion](#conversion) | Unit conversions |
| [sensitivity](#sensitivity) | Sensitivity parameters |
| [noise](#noise) | Noise threshold parameters |
| [iteration](#iteration) | Diffuse iteration parameters |

---

## files

Input file paths for the analysis.

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `datadir` | Path | **Yes** | - | Full path to directory containing FITS files |
| `galaxy` | string | **Yes** | - | Name of dataset (used in output filenames) |
| `starfile` | string | **Yes** | - | Filename of primary stellar/SF tracer FITS file |
| `gasfile` | string | **Yes** | - | Filename of primary gas tracer FITS file |
| `starfile2` | string | No | null | Secondary stellar map for peak ID (if `use_star2=true`) |
| `gasfile2` | string | No | null | Secondary gas map for peak ID (if `use_gas2=true`) |
| `starfile3` | string | No | null | Additional stellar map for flux ratio masking (if `use_star3=true`) |

**Example:**
```yaml
files:
  datadir: /data/galaxies/ngc1234/
  galaxy: NGC1234
  starfile: ngc1234_halpha.fits
  gasfile: ngc1234_co21.fits
```

---

## peak_files

File paths for reusing previously identified peaks (only used when `id_peaks=2`).

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `peaksdir` | Path | No | null | Directory containing .sav files with identified peaks |
| `peakiddir` | Path | No | null | Directory containing clumpfind output files |
| `starpeakidfile` | string | No | null | Peak ID .dat file for stellar map |
| `gaspeakidfile` | string | No | null | Peak ID .dat file for gas map |
| `intpeakidfile` | string | No | null | Interactive peak ID report file |

---

## mask_files

DS9 region files for masking areas in the images.

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `maskdir` | Path | No | null | Directory containing DS9 region files |
| `star_ext_mask` | string | No | null | External mask for starfile (keep regions inside) |
| `star_int_mask` | string | No | null | Internal mask for starfile (exclude regions inside) |
| `gas_ext_mask` | string | No | null | External mask for gasfile |
| `gas_int_mask` | string | No | null | Internal mask for gasfile |
| `star_ext_mask2` | string | No | null | External mask for starfile2 |
| `star_int_mask2` | string | No | null | Internal mask for starfile2 |
| `gas_ext_mask2` | string | No | null | External mask for gasfile2 |
| `gas_int_mask2` | string | No | null | Internal mask for gasfile2 |
| `star_ext_mask3` | string | No | null | External mask for starfile3 |
| `star_int_mask3` | string | No | null | Internal mask for starfile3 |

**Note:** External masks define the region to *include*; internal masks define regions to *exclude*.

---

## unfilt_files

Unfiltered image files for diffuse fraction calculation.

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `unfiltdir` | Path | No | null | Directory containing unfiltered images |
| `star_unfilt_file` | string | No | null | Unfiltered stellar map |
| `gas_unfilt_file` | string | No | null | Unfiltered gas map |

---

## flags1

Module switch flags controlling which analysis steps to run.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mask_images` | bool | true | Apply masks to images |
| `regrid` | bool | true | Read and regrid files, synchronize masks |
| `smoothen` | bool | true | Create smoothed maps for each aperture size |
| `sensitivity` | bool | true | Fit Gaussians to pixel PDFs for sensitivity limits |
| `id_peaks` | 0\|1\|2 | 1 | Peak ID mode: 0=reuse default, 1=identify new, 2=reuse from path |
| `calc_ap_flux` | bool | true | Calculate enclosed flux for each peak/aperture |
| `generate_plot` | bool | true | Generate output plots |
| `get_distances` | bool | true | Calculate distances between peak pairs |
| `calc_obs` | bool | true | Monte-Carlo sample peaks, get observed tuning fork |
| `calc_fit` | bool | true | Fit KL14 uncertainty principle model |
| `diffuse_frac` | bool | true | Calculate diffuse fraction in images |
| `derive_phys` | bool | true | Calculate derived physical quantities |
| `write_output` | bool | true | Write results to output files |
| `cleanup` | 0\|1\|2 | 2 | Cleanup: 0=keep files, 1=prompt, 2=auto-delete |
| `autoexit` | bool | false | Exit automatically on completion |

---

## flags2

Ancillary file usage flags.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_star2` | bool | false | Use different map for stellar peak identification |
| `use_gas2` | bool | false | Use different map for gas peak identification |
| `use_star3` | bool | false | Use additional stellar map for flux ratio masking |

---

## flags3

Masking-related flags.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mstar_ext` | bool | true | Apply external mask to starfile |
| `mstar_int` | bool | true | Apply internal mask to starfile |
| `mgas_ext` | bool | true | Apply external mask to gasfile |
| `mgas_int` | bool | true | Apply internal mask to gasfile |
| `mstar_ext2` | bool | false | Apply external mask to starfile2 |
| `mstar_int2` | bool | false | Apply internal mask to starfile2 |
| `mgas_ext2` | bool | false | Apply external mask to gasfile2 |
| `mgas_int2` | bool | false | Apply internal mask to gasfile2 |
| `mstar_ext3` | bool | false | Apply external mask to starfile3 |
| `mstar_int3` | bool | false | Apply internal mask to starfile3 |
| `convert_masks` | bool | false | Convert masks from other DS9 coordinate systems |
| `cut_radius` | bool | true | Mask outside specified radial interval |

---

## flags4

Analysis option flags.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `set_centre` | bool | true | Use specified centre (true) or image centre (false) |
| `tophat` | bool | true | Use tophat kernel (true) or Gaussian (false) |
| `loglevels` | bool | true | Logarithmic (true) or linear (false) contour spacing |
| `peak_find_tui` | bool | false | Enable interactive peak finding TUI |
| `flux_weight` | bool | false | Flux-weighted position (true) or brightest pixel (false) |
| `calc_ap_area` | bool | true | Calculate aperture area from unmasked pixels |
| `tstar_incl` | bool | false | tstariso includes overlap (true) or tstar=tstariso+tover (false) |
| `peak_prof` | 0\|1\|2 | 2 | Profile: 0=points, 1=constant-density discs, 2=2D Gaussians |
| `map_units` | 0\|1\|2\|3 | 1 | Units: 0=unknown, 1=SFR/gas, 2=gas1/gas2, 3=SFR1/SFR2 |
| `star_tot_mode` | 0\|1 | 0 | Star total: 0=calculate, 1=use star_tot_val |
| `gas_tot_mode` | 0\|1 | 0 | Gas total: 0=calculate, 1=use gas_tot_val |
| `use_X11` | bool | true | Allow X11 windows |
| `log10_output` | bool | true | Write output as log10(value) |

---

## basic_map

Basic map analysis parameters.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `distance` | float | **Required** | Distance to galaxy in **parsecs** |
| `inclination` | float | 0.0 | Inclination angle in degrees (0-90) |
| `posangle` | float | 0.0 | Position angle in degrees |
| `centrex` | int | **Required** | Galaxy centre X pixel (0-indexed from left) |
| `centrey` | int | **Required** | Galaxy centre Y pixel (0-indexed from bottom) |
| `minradius` | float | 0.0 | Minimum radius for analysis in pc |
| `maxradius` | float | **Required** | Maximum radius for analysis in pc |
| `Fs1_Fs2_min` | float | 15.0 | Min flux ratio for pixel inclusion (if use_star3) |
| `max_sample` | int | 10 | Max pixels per resolution FWHM |
| `astr_tolerance` | float | 1.0e-6 | Astrometric tolerance in decimal degrees |
| `nbins` | int | 20 | Bins for sensitivity limit PDF fitting |

**Example:**
```yaml
basic_map:
  distance: 10000000.0   # 10 Mpc in parsecs
  inclination: 30.0
  posangle: 45.0
  centrex: 512
  centrey: 512
  minradius: 0.0
  maxradius: 8000.0
```

---

## aperture

Aperture parameters for the tuning fork analysis.

| Parameter | Type | Default | Constraints | Description |
|-----------|------|---------|-------------|-------------|
| `lapmin` | float | 25.0 | >0 | Minimum aperture diameter in pc |
| `lapmax` | float | 6400.0 | >0 | Maximum aperture diameter in pc |
| `naperture` | int | 9 | >0 | Number of aperture sizes |
| `peak_res` | int | 1 | 0 to naperture-1 | Index for peak ID and min fitting |
| `max_res` | int | 8 | peak_res to naperture-1 | Index for max fitting |

**Notes:**
- Apertures are logarithmically spaced between lapmin and lapmax
- `lapmin` should be at least 2× your resolution FWHM
- `lapmax` should be less than half the field of view

**Example:**
```yaml
aperture:
  lapmin: 50.0
  lapmax: 3200.0
  naperture: 9
  peak_res: 1
  max_res: 7
```

---

## peak_id

Peak identification parameters for the clumpfind algorithm.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `npixmin` | int | 20 | Minimum pixels for a valid peak |
| `nsigma` | float | 5.0 | Significance threshold (× sensitivity limit) |
| `logrange_s` | float | 2.0 | Log range for stellar contour levels |
| `logspacing_s` | float | 0.5 | Log spacing between stellar contours |
| `logrange_g` | float | 2.0 | Log range for gas contour levels |
| `logspacing_g` | float | 0.5 | Log spacing between gas contours |
| `nlinlevel_s` | int | 11 | Linear contour levels for stars (if loglevels=false) |
| `nlinlevel_g` | int | 11 | Linear contour levels for gas (if loglevels=false) |

**Typical values:**
- `npixmin`: 15-50 (lower = more peaks but noisier)
- `nsigma`: 3-10 (lower = more peaks but less significant)

**Example:**
```yaml
peak_id:
  npixmin: 20
  nsigma: 5.0
  logrange_s: 2.0
  logspacing_s: 0.5
  logrange_g: 2.0
  logspacing_g: 0.5
```

---

## timeline

Timescale parameters for the analysis.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tstariso` | float | 1.0 | Reference stellar tracer timescale in Myr |
| `tstariso_errmin` | float | 0.0 | Downward error of tstariso in Myr |
| `tstariso_errmax` | float | 0.0 | Upward error of tstariso in Myr |
| `tgasmini` | float | 0.1 | Minimum tgas during fitting in Myr |
| `tgasmaxi` | float | 1000.0 | Maximum tgas during fitting in Myr |
| `tovermini` | float | 0.01 | Minimum tover during fitting in Myr |

**Typical tstariso values by tracer:**
| Tracer | tstariso (Myr) |
|--------|----------------|
| Hα | 4-5 |
| FUV | 100 |
| 24 μm | 10-20 |
| NUV | 100-200 |

**Example:**
```yaml
timeline:
  tstariso: 4.0          # Hα tracer
  tstariso_errmin: 1.0
  tstariso_errmax: 1.0
  tgasmini: 0.1
  tgasmaxi: 100.0
  tovermini: 0.01
```

---

## fitting

Model fitting parameters.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nmc` | int | 1000 | Monte Carlo peak drawing experiments |
| `ndepth` | int | 4 | Maximum parameter refinement loops |
| `ntry` | int | 101 | Grid size for each parameter |
| `nphysmc` | int | 1000000 | MC experiments for derived physics errors |

**Notes:**
- Increase `nmc` for better error estimates (at cost of runtime)
- `ntry` controls the resolution of the parameter grid
- For quick tests, use `nmc=100`; for production, use `nmc=1000+`

---

## fourier

Fourier filtering parameters.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_unfilt_ims` | bool | false | Use unfiltered images for diffuse fraction |
| `diffuse_quant` | 0\|1 | 1 | Diffuse quantity: 0=flux, 1=power |
| `f_filter_type` | 0\|1\|2 | 2 | Filter: 0=butterworth, 1=gaussian, 2=ideal |
| `bw_order` | int | 2 | Butterworth filter order |
| `filter_len_conv` | float | 1.0 | Filter cut length = λ × filter_len_conv |
| `emfrac_cor_mode` | 0-5 | 0 | Correction: 0=none, 1=flux loss, 4=overlap, 5=both |
| `rpeak_cor_mode` | 0\|1 | 0 | rpeak: 0=measured, 1=supplied values |
| `rpeaks_cor_val` | float | 1.0 | r_peak_star for flux loss correction |
| `rpeaks_cor_emin` | float | 0.5 | Downward error of r_peak_star |
| `rpeaks_cor_emax` | float | 0.5 | Upward error of r_peak_star |
| `rpeakg_cor_val` | float | 1.0 | r_peak_gas for flux loss correction |
| `rpeakg_cor_emin` | float | 0.5 | Downward error of r_peak_gas |
| `rpeakg_cor_emax` | float | 0.5 | Upward error of r_peak_gas |

---

## conversion

Unit conversion parameters.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `convstar` | float | -3.69206 | log10 conversion: pixel → SFR (M☉/yr) or mass (M☉) |
| `convstar_rerr` | float | 0.0 | Relative error of convstar |
| `convgas` | float | 2.30794 | log10 conversion: pixel → mass (M☉) or SFR (M☉/yr) |
| `convgas_rerr` | float | 0.0 | Relative error of convgas |
| `convstar3` | float | 0.0 | log10 conversion for starfile3 |
| `convstar3_rerr` | float | 0.0 | Relative error of convstar3 |
| `lighttomass` | float | 0.002 | Light-to-mass ratio in m² s⁻³ (default: SNe) |
| `momratetomass` | float | 5.0e-10 | Momentum rate per mass in m s⁻² |
| `star_tot_val` | float | 1.0 | Total SFR/mass for stellar map (if star_tot_mode=1) |
| `star_tot_err` | float | 0.1 | Error on star_tot_val |
| `gas_tot_val` | float | 1.0e9 | Total mass/SFR for gas map (if gas_tot_mode=1) |
| `gas_tot_err` | float | 1.0e8 | Error on gas_tot_val |

---

## sensitivity

Sensitivity parameters.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_stds` | bool | false | Use supplied std values (true) or calculate (false) |
| `std_star` | float | 0.1 | Standard deviation of starfile |
| `std_star3` | float | 0.1 | Standard deviation of starfile3 |
| `std_gas` | float | 0.1 | Standard deviation of gasfile |

---

## noise

Noise threshold parameters.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_noisecut` | bool | true | Mask below noise threshold after filtering |
| `noisethresh_s` | float | 20.0 | Noise threshold for star map |
| `noisethresh_g` | float | 20.0 | Noise threshold for gas map |

---

## iteration

Iterative diffuse filtering parameters.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_guess` | bool | false | Pre-filter with initial λ guess |
| `initial_guess` | float | 200.0 | Initial λ estimate in pc |
| `iter_criterion` | float | 0.05 | Fractional convergence criterion |
| `iter_crit_len` | int | 2 | Previous iterations for convergence check |
| `iter_nmax` | int | 10 | Maximum iterations |
| `iter_filter` | 0\|1\|2 | 2 | Filter: 0=butterworth, 1=gaussian, 2=ideal |
| `iter_bwo` | int | 2 | Butterworth order for iteration |
| `iter_len_conv` | float | 2.0 | Filter cut length conversion factor |
| `iter_rpeak_mode` | 0\|1 | 0 | rpeak: 0=use rpeak_cor_mode, 1=use iter0 |
| `iter_tot_mode_s` | 0\|1\|2 | 0 | Star total: 0=per-iter, 1=iter0, 2=star_tot_val |
| `iter_tot_mode_g` | 0\|1\|2 | 0 | Gas total: 0=per-iter, 1=iter0, 2=gas_tot_val |
| `iter_autoexit` | bool | false | Auto-exit on convergence |
| `use_nice` | bool | false | Use nice for process priority |
| `nice_value` | int | 0 | Nice value (-20 to 19) |

**Example:**
```yaml
iteration:
  use_guess: false
  initial_guess: 200.0
  iter_criterion: 0.05
  iter_nmax: 10
  iter_filter: 2
```

---

## Complete Example Configuration

```yaml
# Heisenberg Configuration File
# Galaxy: NGC1234

files:
  datadir: /data/galaxies/ngc1234/
  galaxy: NGC1234
  starfile: ngc1234_halpha.fits
  gasfile: ngc1234_co21.fits

mask_files:
  maskdir: /data/galaxies/ngc1234/masks/
  star_ext_mask: ngc1234_field.reg
  gas_ext_mask: ngc1234_field.reg

flags1:
  mask_images: true
  id_peaks: 1
  calc_fit: true
  diffuse_frac: true

flags4:
  tophat: true
  loglevels: true
  peak_prof: 2

basic_map:
  distance: 10000000.0
  inclination: 30.0
  posangle: 45.0
  centrex: 512
  centrey: 512
  maxradius: 8000.0

aperture:
  lapmin: 50.0
  lapmax: 3200.0
  naperture: 9

peak_id:
  npixmin: 20
  nsigma: 5.0
  logrange_s: 2.0
  logspacing_s: 0.5
  logrange_g: 2.0
  logspacing_g: 0.5

timeline:
  tstariso: 4.0
  tstariso_errmin: 1.0
  tstariso_errmax: 1.0
  tgasmini: 0.1
  tgasmaxi: 100.0
  tovermini: 0.01

fitting:
  nmc: 1000
  ndepth: 4
  ntry: 101

iteration:
  iter_criterion: 0.05
  iter_nmax: 10
```

---

## See Also

- [Quick Start](../quickstart.md) - Getting started guide
- [CLI Reference](cli-reference.md) - Command-line usage
- [Parameter Tuning](../tutorials/parameter-tuning.md) - Best practices for parameter selection

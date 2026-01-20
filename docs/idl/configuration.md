# IDL Configuration Format

Complete specification of the IDL input file format.

## Overview

The IDL input file is a plain text file with parameters listed in a specific order. Each parameter occupies one line, with optional comments after semicolons.

**File structure:**
```
#==== SECTION HEADER ====#
value1                   ; parameter1 description
value2                   ; parameter2 description
...
```

## Complete Input File Specification

### Section 1: File Names

```
#==== INPUT FILE NAME ====#
/full/path/to/data/      ; datadir - directory containing FITS files
NGC1234                  ; galaxy - dataset name for output files
star.fits                ; starfile - stellar/SF tracer FITS file
star2.fits               ; starfile2 - secondary stellar map (or 0 if unused)
gas.fits                 ; gasfile - gas tracer FITS file
gas2.fits                ; gasfile2 - secondary gas map (or 0 if unused)
star3.fits               ; starfile3 - additional stellar map (or 0 if unused)
```

### Section 2: Peak ID File Names

```
#==== PEAK ID FILE NAME ====#
/path/to/peaks/          ; peaksdir - directory with .sav peak files
/path/to/peakid/         ; peakiddir - directory with clumpfind output
star_peaks.dat           ; starpeakidfile - stellar peak ID file
gas_peaks.dat            ; gaspeakidfile - gas peak ID file
interactive_report.dat   ; intpeakidfile - interactive peak report
```

### Section 3: Mask File Names

```
#==== MASK FILE NAME ====#
/path/to/masks/          ; maskdir - directory with DS9 region files
star_ext.reg             ; star_ext_mask - external mask for starfile
star_int.reg             ; star_int_mask - internal mask for starfile
gas_ext.reg              ; gas_ext_mask - external mask for gasfile
gas_int.reg              ; gas_int_mask - internal mask for gasfile
star2_ext.reg            ; star_ext_mask2
star2_int.reg            ; star_int_mask2
gas2_ext.reg             ; gas_ext_mask2
gas2_int.reg             ; gas_int_mask2
star3_ext.reg            ; star_ext_mask3
star3_int.reg            ; star_int_mask3
```

### Section 4: Unfiltered File Names

```
#==== UNFILTERED FILE NAME ====#
/path/to/unfiltered/     ; unfiltdir
star_unfilt.fits         ; star_unfilt_file
gas_unfilt.fits          ; gas_unfilt_file
```

### Section 5: Flags 1 (Module Switches)

```
#==== FLAGS 1 ====#
1 1 1 1 1 1 1 1 1 1 1 1 1 2 0
; mask_images regrid smoothen sensitivity id_peaks calc_ap_flux
; generate_plot get_distances calc_obs calc_fit diffuse_frac
; derive_phys write_output cleanup autoexit
```

| Position | Parameter | Values | Description |
|----------|-----------|--------|-------------|
| 1 | mask_images | 0/1 | Apply masks |
| 2 | regrid | 0/1 | Read and regrid files |
| 3 | smoothen | 0/1 | Create smoothed maps |
| 4 | sensitivity | 0/1 | Calculate sensitivity |
| 5 | id_peaks | 0/1/2 | Peak ID mode |
| 6 | calc_ap_flux | 0/1 | Calculate aperture flux |
| 7 | generate_plot | 0/1 | Generate plots |
| 8 | get_distances | 0/1 | Calculate distances |
| 9 | calc_obs | 0/1 | Monte Carlo sampling |
| 10 | calc_fit | 0/1 | Fit KL14 model |
| 11 | diffuse_frac | 0/1 | Calculate diffuse fraction |
| 12 | derive_phys | 0/1 | Derive physics |
| 13 | write_output | 0/1 | Write output |
| 14 | cleanup | 0/1/2 | Cleanup mode |
| 15 | autoexit | 0/1 | Auto-exit |

### Section 6: Flags 2 (Ancillary Files)

```
#==== FLAGS 2 ====#
0 0 0
; use_star2 use_gas2 use_star3
```

### Section 7: Flags 3 (Masking)

```
#==== FLAGS 3 ====#
1 1 1 1 0 0 0 0 0 0 0 1
; mstar_ext mstar_int mgas_ext mgas_int
; mstar_ext2 mstar_int2 mgas_ext2 mgas_int2
; mstar_ext3 mstar_int3 convert_masks cut_radius
```

### Section 8: Flags 4 (Analysis Options)

```
#==== FLAGS 4 ====#
1 1 1 0 0 1 0 2 1 0 0 1 1
; set_centre tophat loglevels peak_find_tui flux_weight
; calc_ap_area tstar_incl peak_prof map_units
; star_tot_mode gas_tot_mode use_X11 log10_output
```

| Position | Parameter | Values | Description |
|----------|-----------|--------|-------------|
| 1 | set_centre | 0/1 | Use specified centre |
| 2 | tophat | 0/1 | Tophat (1) or Gaussian (0) kernel |
| 3 | loglevels | 0/1 | Logarithmic contours |
| 4 | peak_find_tui | 0/1 | Interactive peak finding |
| 5 | flux_weight | 0/1 | Flux-weighted positions |
| 6 | calc_ap_area | 0/1 | Calculate aperture area |
| 7 | tstar_incl | 0/1 | tstar includes overlap |
| 8 | peak_prof | 0/1/2 | Peak profile model |
| 9 | map_units | 0/1/2/3 | Map unit type |
| 10 | star_tot_mode | 0/1 | Star total mode |
| 11 | gas_tot_mode | 0/1 | Gas total mode |
| 12 | use_X11 | 0/1 | Allow X11 windows |
| 13 | log10_output | 0/1 | Log10 output values |

### Section 9: Basic Map Parameters

```
#==== INPUT PARAMETERS 1 ====#
10000000.0               ; distance [pc]
30.0                     ; inclination [degrees]
45.0                     ; posangle [degrees]
512                      ; centrex [pixels]
512                      ; centrey [pixels]
0.0                      ; minradius [pc]
8000.0                   ; maxradius [pc]
15.0                     ; Fs1_Fs2_min
10                       ; max_sample
1.0e-6                   ; astr_tolerance [degrees]
20                       ; nbins
```

### Section 10: Aperture Parameters

```
#==== INPUT PARAMETERS 2 ====#
50.0                     ; lapmin [pc]
3200.0                   ; lapmax [pc]
9                        ; naperture
1                        ; peak_res (0-indexed)
7                        ; max_res (0-indexed)
```

### Section 11: Peak ID Parameters

```
#==== INPUT PARAMETERS 3 ====#
20                       ; npixmin
5.0                      ; nsigma
2.0                      ; logrange_s
0.5                      ; logspacing_s
2.0                      ; logrange_g
0.5                      ; logspacing_g
11                       ; nlinlevel_s
11                       ; nlinlevel_g
```

### Section 12: Timeline Parameters

```
#==== INPUT PARAMETERS 4 ====#
4.0                      ; tstariso [Myr]
1.0                      ; tstariso_errmin [Myr]
1.0                      ; tstariso_errmax [Myr]
0.1                      ; tgasmini [Myr]
100.0                    ; tgasmaxi [Myr]
0.01                     ; tovermini [Myr]
```

### Section 13: Fitting Parameters

```
#==== INPUT PARAMETERS 5 ====#
1000                     ; nmc
4                        ; ndepth
101                      ; ntry
1000000                  ; nphysmc
```

### Section 14: Fourier Filter Parameters

```
#==== INPUT PARAMETERS 6 ====#
0                        ; use_unfilt_ims
1                        ; diffuse_quant (0=flux, 1=power)
2                        ; f_filter_type (0=bw, 1=gauss, 2=ideal)
2                        ; bw_order
1.0                      ; filter_len_conv
0                        ; emfrac_cor_mode
0                        ; rpeak_cor_mode
1.0                      ; rpeaks_cor_val
0.5                      ; rpeaks_cor_emin
0.5                      ; rpeaks_cor_emax
1.0                      ; rpeakg_cor_val
0.5                      ; rpeakg_cor_emin
0.5                      ; rpeakg_cor_emax
```

### Section 15: Conversion Parameters

```
#==== INPUT PARAMETERS 7 ====#
-3.69206                 ; convstar (log10)
0.0                      ; convstar_rerr
2.30794                  ; convgas (log10)
0.0                      ; convgas_rerr
0.0                      ; convstar3
0.0                      ; convstar3_rerr
0.002                    ; lighttomass [m^2 s^-3]
5.0e-10                  ; momratetomass [m s^-2]
1.0                      ; star_tot_val
0.1                      ; star_tot_err
1.0e9                    ; gas_tot_val
1.0e8                    ; gas_tot_err
```

### Section 16: Sensitivity Parameters

```
#==== INPUT PARAMETERS 8 ====#
0                        ; use_stds
0.1                      ; std_star
0.1                      ; std_star3
0.1                      ; std_gas
```

### Section 17: Noise Threshold Parameters

```
#==== INPUT PARAMETERS 9 ====#
1                        ; use_noisecut
20.0                     ; noisethresh_s
20.0                     ; noisethresh_g
```

### Section 18: Iteration Parameters

```
#==== INPUT PARAMETERS 10 ====#
0                        ; use_guess
200.0                    ; initial_guess [pc]
0.05                     ; iter_criterion
2                        ; iter_crit_len
10                       ; iter_nmax
2                        ; iter_filter (0=bw, 1=gauss, 2=ideal)
2                        ; iter_bwo
2.0                      ; iter_len_conv
0                        ; iter_rpeak_mode
0                        ; iter_tot_mode_s
0                        ; iter_tot_mode_g
0                        ; iter_autoexit
0                        ; use_nice
0                        ; nice_value
```

## Complete Example Input File

```
#==== INPUT FILE NAME ====#
/data/galaxies/ngc1234/
NGC1234
ngc1234_halpha.fits
0
ngc1234_co21.fits
0
0
#==== PEAK ID FILE NAME ====#
0
0
0
0
0
#==== MASK FILE NAME ====#
/data/galaxies/ngc1234/masks/
ngc1234_field.reg
0
ngc1234_field.reg
0
0
0
0
0
0
0
#==== UNFILTERED FILE NAME ====#
0
0
0
#==== FLAGS 1 ====#
1 1 1 1 1 1 1 1 1 1 1 1 1 2 0
#==== FLAGS 2 ====#
0 0 0
#==== FLAGS 3 ====#
1 0 1 0 0 0 0 0 0 0 0 1
#==== FLAGS 4 ====#
1 1 1 0 0 1 0 2 1 0 0 1 1
#==== INPUT PARAMETERS 1 ====#
10000000.0
30.0
45.0
512
512
0.0
8000.0
15.0
10
1.0e-6
20
#==== INPUT PARAMETERS 2 ====#
50.0
3200.0
9
1
7
#==== INPUT PARAMETERS 3 ====#
20
5.0
2.0
0.5
2.0
0.5
11
11
#==== INPUT PARAMETERS 4 ====#
4.0
1.0
1.0
0.1
100.0
0.01
#==== INPUT PARAMETERS 5 ====#
1000
4
101
1000000
#==== INPUT PARAMETERS 6 ====#
0
1
2
2
1.0
0
0
1.0
0.5
0.5
1.0
0.5
0.5
#==== INPUT PARAMETERS 7 ====#
-3.69206
0.0
2.30794
0.0
0.0
0.0
0.002
5.0e-10
1.0
0.1
1.0e9
1.0e8
#==== INPUT PARAMETERS 8 ====#
0
0.1
0.1
0.1
#==== INPUT PARAMETERS 9 ====#
1
20.0
20.0
#==== INPUT PARAMETERS 10 ====#
0
200.0
0.05
2
10
2
2
2.0
0
0
0
0
0
0
```

## Notes

### Using 0 for Unused Parameters

For optional file paths, use `0` to indicate the parameter is not used:
```
0                        ; starfile2 - not using secondary stellar map
```

### Comments

Comments can be added after values using semicolons:
```
10000000.0               ; distance in pc (10 Mpc)
```

### Paths

- `datadir` must be an absolute path ending with `/`
- Other paths are relative to their respective directories
- Use `0` for unused directory parameters

### Boolean Values

- Use `1` for true/yes/enabled
- Use `0` for false/no/disabled

### Creating Input Files

Use the IDL helper to create a template:
```idl
IDL> make_input_file, '/path/to/new_input_file'
```

## See Also

- [IDL Usage Guide](usage.md) - Running the analysis
- [Configuration Reference](../python/configuration.md) - Parameter descriptions
- [Parameter Tuning](../tutorials/parameter-tuning.md) - Best practices

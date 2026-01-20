# IDL Usage Guide

Complete guide to running Heisenberg using the IDL implementation.

## Overview

The IDL implementation provides two main entry points:

| Entry Point | Description |
|-------------|-------------|
| `heisenberg.pro` | Full analysis with iterative diffuse filtering (recommended) |
| `heisenberg_nodf.pro` | Single-pass analysis without diffuse filtering |

## Prerequisites

- **IDL 8.0+** (Harris Geospatial Solutions)
- **IDL Astronomy User's Library** ([download](https://idlastro.gsfc.nasa.gov/))
- Input FITS files with matching WCS

## Quick Start

### Running with Diffuse Filtering (Recommended)

```idl
; From command line
idl heisenberg -arg /path/to/input_file

; Or interactively in IDL
IDL> inputfile = '/path/to/input_file'
IDL> .run heisenberg
```

### Running without Diffuse Filtering

```idl
; From command line
idl heisenberg_nodf -arg /path/to/input_file

; Or interactively
IDL> inputfile = '/path/to/input_file'
IDL> .run heisenberg_nodf
```

## Execution Modes

### 1. Iterative Diffuse Filtering (`heisenberg.pro`)

This is the recommended mode for scientific analyses. The procedure:

1. Loads configuration and FITS data
2. Identifies peaks in stellar and gas maps
3. Runs initial tuningfork fit to get λ
4. Filters diffuse emission at scale λ using Fourier methods
5. Re-identifies peaks and refits
6. Repeats until λ converges (within `iter_criterion`)

**Key iteration parameters:**
- `iter_criterion`: Fractional convergence threshold (default: 0.05)
- `iter_nmax`: Maximum iterations (default: 10)
- `iter_filter`: Filter type (0=Butterworth, 1=Gaussian, 2=Ideal)

### 2. Single-Pass Analysis (`heisenberg_nodf.pro`)

Runs a single tuningfork fit without diffuse filtering. Useful for:
- Quick parameter exploration
- Pre-filtered data
- Comparing results with/without filtering

### 3. Interactive Peak Finding

Enable the interactive TUI for peak finding by setting in your input file:

```
1                    ; peak_find_tui (FLAGS 4, line 4)
```

This allows you to:
- Visualize detected peaks overlaid on maps
- Adjust `npixmin` and `nsigma` interactively
- Modify contour levels
- Accept or reject peak configurations

## Input File Format

The IDL input file is a text file with parameters in a specific order. See [IDL Configuration](configuration.md) for the complete format specification.

**Example minimal input file structure:**
```
#==== INPUT FILE NAME ====#
/path/to/data/           ; datadir
NGC1234                  ; galaxy name
star.fits                ; starfile
...                      ; (continue with all parameters)
```

## Output Files

The analysis produces several output files in the data directory:

| File Pattern | Description |
|--------------|-------------|
| `<galaxy>_output.txt` | Main results summary |
| `<galaxy>_output_log.txt` | Detailed log |
| `<galaxy>_tuningfork.ps` | Tuning fork diagram (PostScript) |
| `<galaxy>_star_peaks.dat` | Stellar peak positions |
| `<galaxy>_gas_peaks.dat` | Gas peak positions |
| `<galaxy>_star_peaks.reg` | DS9 region file for stellar peaks |
| `<galaxy>_gas_peaks.reg` | DS9 region file for gas peaks |

### Reading Output in IDL

```idl
; Read the output file
output = read_heisenberg_output('/path/to/NGC1234_output.txt')

; Access results
print, 'tgas = ', output.tgas, ' +', output.tgas_errmax, ' -', output.tgas_errmin
print, 'tover = ', output.tover, ' +', output.tover_errmax, ' -', output.tover_errmin
print, 'lambda = ', output.lambda, ' +', output.lambda_errmax, ' -', output.lambda_errmin
```

## Command-Line Execution

### Basic Usage

```bash
# With diffuse filtering
idl heisenberg -arg /absolute/path/to/input_file

# Without diffuse filtering
idl heisenberg_nodf -arg /absolute/path/to/input_file
```

### Batch Processing

```bash
#!/bin/bash
# Process multiple galaxies

for input in /data/configs/*.inp; do
    echo "Processing: $input"
    idl heisenberg -arg "$input"
done
```

### Background Execution

```bash
# Run in background with output logging
nohup idl heisenberg -arg /path/to/input_file > heisenberg.log 2>&1 &
```

## Interactive Session

For exploratory analysis, work interactively in IDL:

```idl
; Start IDL
IDL>

; Set input file path
IDL> inputfile = '/path/to/input_file'

; Run the analysis
IDL> .run heisenberg

; After completion, explore results
IDL> help, /structures  ; See available structures
IDL> print, tgas, tover, lambda
```

## Module Flags

Control which analysis steps run via FLAGS 1 in the input file:

| Flag | Name | Description |
|------|------|-------------|
| 1 | mask_images | Apply DS9 region masks |
| 2 | regrid | Read and regrid FITS files |
| 3 | smoothen | Create smoothed aperture maps |
| 4 | sensitivity | Calculate sensitivity limits |
| 5 | id_peaks | Identify peaks (0/1/2) |
| 6 | calc_ap_flux | Calculate aperture fluxes |
| 7 | generate_plot | Generate PostScript plots |
| 8 | get_distances | Calculate peak pair distances |
| 9 | calc_obs | Monte Carlo sampling |
| 10 | calc_fit | Fit KL14 model |
| 11 | diffuse_frac | Calculate diffuse fraction |
| 12 | derive_phys | Derive physical quantities |
| 13 | write_output | Write output files |
| 14 | cleanup | Clean temporary files (0/1/2) |
| 15 | autoexit | Exit on completion |

**Example:** To skip plotting and run non-interactively:
```
1 1 1 1 1 1 0 1 1 1 1 1 1 2 1    ; FLAGS 1
```

## Visualization with DS9

After running, visualize results in DS9:

```bash
# Open stellar map with peak regions
ds9 /path/to/star.fits -region /path/to/NGC1234_star_peaks.reg

# Or load both maps
ds9 -multiframe \
    /path/to/star.fits -region NGC1234_star_peaks.reg \
    /path/to/gas.fits -region NGC1234_gas_peaks.reg
```

## Common Workflows

### First Analysis of New Galaxy

1. Create input file from template:
```idl
IDL> make_input_file, '/path/to/new_input_file'
```

2. Edit input file with your parameters

3. Run with interactive peak finding to tune parameters:
```idl
IDL> inputfile = '/path/to/new_input_file'
IDL> .run heisenberg
; When prompted, adjust peak finding parameters
```

4. Review output and iterate

### Re-running with Different Parameters

To reuse previously identified peaks:

1. Set `id_peaks = 0` in FLAGS 1 (reuse from default location)
2. Or set `id_peaks = 2` and specify peak file paths

```
0 1 1 1 0 1 1 1 1 1 1 1 1 2 0    ; FLAGS 1 (id_peaks=0)
```

### Comparing With/Without Diffuse Filtering

```bash
# Run without filtering first
idl heisenberg_nodf -arg input_file
mv NGC1234_output.txt NGC1234_output_nodf.txt

# Run with filtering
idl heisenberg -arg input_file

# Compare results
diff NGC1234_output_nodf.txt NGC1234_output.txt
```

## Troubleshooting

### "Procedure not found"

Ensure IDL path includes the Heisenberg directory:
```idl
IDL> !PATH = !PATH + ':' + expand_path('+/path/to/heisenberg/IDL/')
```

### "FITS file not found"

- Check `datadir` path is absolute and correct
- Verify file permissions
- Ensure filenames match exactly (case-sensitive)

### "Astrometry mismatch"

Star and gas maps must have matching WCS:
```idl
; Check WCS in IDL
IDL> star_header = headfits('/path/to/star.fits')
IDL> gas_header = headfits('/path/to/gas.fits')
IDL> print, sxpar(star_header, 'CRVAL1'), sxpar(gas_header, 'CRVAL1')
```

### Memory Issues

For large images, reduce `max_sample` parameter or process a smaller region using `minradius`/`maxradius`.

### X11 Display Errors

If running remotely without X11:
```
0                    ; use_X11 = 0 in FLAGS 4
```

Or set up X11 forwarding:
```bash
ssh -X user@remote
```

## See Also

- [IDL Configuration](configuration.md) - Input file format specification
- [IDL Interactive Mode](interactive-mode.md) - Using the peak finding TUI
- [IDL Procedures Reference](procedures.md) - All IDL procedures
- [Installation Guide](../installation.md) - IDL setup instructions

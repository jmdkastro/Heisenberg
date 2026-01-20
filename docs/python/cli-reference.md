# CLI Reference

Complete reference for the Heisenberg command-line interface.

## Overview

Heisenberg provides a command-line interface via the `heisenberg` command. The CLI is built using [Click](https://click.palletsprojects.com/) and supports multiple subcommands.

```bash
heisenberg [OPTIONS] COMMAND [ARGS]...
```

## Global Options

| Option | Description |
|--------|-------------|
| `--version` | Show version number and exit |
| `--help` | Show help message and exit |

## Commands

### heisenberg run

Run the Heisenberg analysis on observational data.

```bash
heisenberg run INPUT_FILE [OPTIONS]
```

#### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `INPUT_FILE` | Yes | Path to configuration file (YAML or IDL input_file format) |

#### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--no-diffuse-filtering` | `-n` | Flag | False | Skip iterative diffuse filtering; run single-pass analysis |
| `--interactive` | `-i` | Flag | False | Launch interactive peak finding TUI before analysis |
| `--output` | `-o` | Path | `<galaxy>_results.json` | Output file path for results |
| `--output-dir` | `-d` | Path | `.` | Output directory for plots and region files |
| `--plot` | `-p` | Flag | False | Generate diagnostic plots (tuning fork, PDFs) |
| `--regions` | `-r` | Flag | False | Export peak positions to DS9 region files |
| `--verbose` | `-v` | Flag | False | Enable verbose output with progress information |

#### Examples

**Basic analysis with diffuse filtering (recommended):**
```bash
heisenberg run my_config.yaml
```

**Single-pass analysis without diffuse filtering:**
```bash
heisenberg run my_config.yaml --no-diffuse-filtering
```

**Full analysis with all outputs:**
```bash
heisenberg run my_config.yaml --plot --regions --verbose --output results/ngc1234.json --output-dir results/
```

**Interactive peak finding:**
```bash
heisenberg run my_config.yaml --interactive --verbose
```

**Using short options:**
```bash
heisenberg run my_config.yaml -n -p -v -o output.json
```

#### Output Files

When run completes, the following files are created:

| File | Created When | Description |
|------|--------------|-------------|
| `<galaxy>_results.json` | Always | Full results in JSON format |
| `<galaxy>_results.txt` | Always | Human-readable summary |
| `<galaxy>_tuningfork.png` | `--plot` | Tuning fork diagram |
| `<galaxy>_star_peaks.reg` | `--regions` | Stellar peak positions (DS9 format) |
| `<galaxy>_gas_peaks.reg` | `--regions` | Gas peak positions (DS9 format) |

#### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (configuration, data loading, or analysis failure) |

---

### heisenberg validate

Validate a configuration file without running the analysis.

```bash
heisenberg validate INPUT_FILE
```

#### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `INPUT_FILE` | Yes | Path to configuration file to validate |

#### Examples

```bash
heisenberg validate my_config.yaml
```

**Successful output:**
```
Configuration is valid!
  Galaxy: NGC1234
  Distance: 10000000 pc
  Data directory: /path/to/data/
  Star file: ngc1234_halpha.fits
  Gas file: ngc1234_co.fits
  Apertures: 9 (50 - 3200 pc)
```

**Error output:**
```
Configuration error: Missing required field 'files.starfile'
```

#### What's Validated

- All required parameters are present
- Parameter types are correct (numbers, booleans, paths)
- File paths exist (if datadir is accessible)
- Parameter values are within valid ranges
- Logical consistency (e.g., lapmin < lapmax)

---

### heisenberg create-input

Create a template configuration file.

```bash
heisenberg create-input OUTPUT_FILE [OPTIONS]
```

#### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `OUTPUT_FILE` | Yes | Path where template will be written |

#### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--template` | `-t` | Choice | `full` | Template type: `minimal` (required only) or `full` (all parameters) |

#### Examples

**Create full template:**
```bash
heisenberg create-input my_config.yaml
```

**Create minimal template:**
```bash
heisenberg create-input my_config.yaml --template minimal
```

---

### heisenberg quick-test

Run a quick test with synthetic data to verify installation.

```bash
heisenberg quick-test [OPTIONS]
```

#### Options

| Option | Short | Type | Default | Description |
|--------|-------|------|---------|-------------|
| `--n-peaks` | `-n` | Integer | 15 | Number of synthetic peaks to generate |
| `--output` | `-o` | Path | None | Save results to file |
| `--plot` | `-p` | Flag | False | Generate and display diagnostic plot |

#### Examples

**Basic test:**
```bash
heisenberg quick-test
```

**Test with plot:**
```bash
heisenberg quick-test --plot
```

**Test with more peaks and save results:**
```bash
heisenberg quick-test --n-peaks 30 --output test_results.json --plot
```

#### What It Does

1. Generates synthetic star and gas images with Gaussian peaks
2. Runs the tuningfork analysis pipeline
3. Displays fitted parameters (t_gas, t_over, lambda)
4. Optionally saves results and generates plots

This is useful for:
- Verifying installation is working
- Understanding output format
- Testing without real data

---

## Configuration File Formats

The CLI accepts two configuration file formats:

### YAML Format (Recommended)

```yaml
files:
  datadir: /path/to/data/
  galaxy: NGC1234
  starfile: star.fits
  gasfile: gas.fits

basic_map:
  distance: 10000000.0
  inclination: 30.0

aperture:
  lapmin: 50.0
  lapmax: 3200.0
  naperture: 9
```

See [Configuration Reference](configuration.md) for all parameters.

### IDL Input File Format

```
#==== INPUT FILE NAME ====#
/path/to/data/
NGC1234
star.fits
gas.fits
...
```

This format matches the original IDL code and is supported for backward compatibility.

---

## Environment Variables

| Variable | Description |
|----------|-------------|
| `HEISENBERG_CONFIG` | Default configuration file path |
| `HEISENBERG_OUTPUT_DIR` | Default output directory |

---

## Shell Completion

Enable shell completion for bash:

```bash
# Add to ~/.bashrc
eval "$(_HEISENBERG_COMPLETE=bash_source heisenberg)"
```

For zsh:

```bash
# Add to ~/.zshrc
eval "$(_HEISENBERG_COMPLETE=zsh_source heisenberg)"
```

---

## Common Workflows

### First-time analysis of a new galaxy

```bash
# 1. Create configuration from template
heisenberg create-input ngc1234_config.yaml

# 2. Edit configuration with your data paths and parameters
# (use your favorite editor)

# 3. Validate configuration
heisenberg validate ngc1234_config.yaml

# 4. Run with interactive peak finding to tune parameters
heisenberg run ngc1234_config.yaml --interactive --verbose

# 5. Run full analysis
heisenberg run ngc1234_config.yaml --plot --regions --verbose
```

### Quick exploration / testing

```bash
# Single-pass analysis without diffuse filtering
heisenberg run config.yaml -n -v
```

### Production analysis

```bash
# Full analysis with all outputs
heisenberg run config.yaml \
    --output results/ngc1234_results.json \
    --output-dir results/ \
    --plot \
    --regions \
    --verbose
```

### Batch processing

```bash
# Process multiple galaxies
for config in configs/*.yaml; do
    galaxy=$(basename "$config" .yaml)
    heisenberg run "$config" \
        --output "results/${galaxy}.json" \
        --output-dir "results/${galaxy}/" \
        --plot
done
```

---

## Troubleshooting

### Command not found

Ensure the package is installed and the Python scripts directory is in your PATH:

```bash
pip install -e Python/
which heisenberg
```

### Permission denied

Check file permissions on input data and output directory:

```bash
ls -la /path/to/data/
mkdir -p output_dir && chmod 755 output_dir
```

### Memory errors

For large images, reduce `max_sample` or process in chunks:

```yaml
basic_map:
  max_sample: 5  # Reduce from default 10
```

### See Also

- [Configuration Reference](configuration.md) - All configuration parameters
- [Troubleshooting](../troubleshooting.md) - Common issues and solutions
- [Quick Start](../quickstart.md) - Getting started guide

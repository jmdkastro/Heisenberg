# Glossary

This glossary defines key terms and symbols used in Heisenberg and the KL14 uncertainty principle methodology.

## Core Parameters

| Term | Symbol | Units | Definition |
|------|--------|-------|------------|
| **Gas lifetime** | t_gas | Myr | The duration from the onset of gas cloud assembly to complete dispersal. This is the total time a molecular cloud exists as a coherent structure. |
| **Overlap time** | t_over | Myr | The duration during which both the gas and stellar tracers are simultaneously visible. This corresponds to the embedded phase of star formation before feedback disperses the gas. |
| **Stellar timescale** | t_star | Myr | The visibility timescale of the stellar tracer. This is an **input parameter** (not fitted) that depends on the tracer used (e.g., ~4-5 Myr for Hα, ~100 Myr for FUV). |
| **Total timeline** | t_total | Myr | The total evolutionary timeline: t_total = t_gas + t_star - t_over |
| **Region separation length** | λ (lambda) | pc | The characteristic spacing between independent star-forming regions. Regions separated by more than λ evolve independently and are decorrelated in time. |

## Derived Quantities

| Term | Symbol | Units | Definition |
|------|--------|-------|------------|
| **Star formation efficiency** | ε_sf | dimensionless | The fraction of gas converted to stars per cloud lifecycle. Calculated from the fitted timescales. |
| **Feedback timescale** | t_fb | Myr | Time for feedback to disrupt the parent cloud after star formation begins. Related to t_over. |
| **Feedback velocity** | v_fb | km/s | Characteristic velocity at which feedback-driven material expands. |
| **Feedback outflow rate** | M_dot_fb | M☉/yr | Mass outflow rate driven by stellar feedback. |

## Analysis Concepts

### Tuning Fork Diagram

The characteristic plot showing the gas-to-stellar flux ratio (or stellar-to-gas flux ratio) as a function of aperture size. Named for its fork-like shape:

- **At small apertures**: Large scatter as individual regions are resolved; flux ratios depend on the evolutionary phase of each region
- **At large apertures**: Ratios converge to the galactic average as many regions are averaged together
- **Branch separation**: The two branches (gas-dominated and star-dominated regions) merge at scales larger than λ

### Diffuse Emission

Extended emission that is not associated with discrete peaks/clumps. Sources include:

- Old stellar populations (for stellar tracers)
- Unresolved distant sources
- Truly diffuse gas not associated with star-forming regions
- Instrumental effects (PSF wings, scattered light)

Diffuse emission is removed via Fourier filtering at the scale λ.

### Peak / Clump

A localized overdensity in the stellar or gas map identified by the clumpfind algorithm. Peaks represent individual star-forming regions or GMCs.

## Configuration Terms

### Files Group

| Term | Description |
|------|-------------|
| **datadir** | Directory containing the input FITS files |
| **starfile** | FITS file containing the stellar/star formation tracer map |
| **gasfile** | FITS file containing the gas tracer map |
| **galaxy** | Identifier string used in output filenames |

### Aperture Parameters

| Term | Description |
|------|-------------|
| **lapmin** | Minimum aperture diameter in parsecs |
| **lapmax** | Maximum aperture diameter in parsecs |
| **naperture** | Number of aperture sizes (logarithmically spaced between lapmin and lapmax) |
| **peak_res** | Index of aperture size used for peak identification (0-indexed) |

### Peak Detection Parameters

| Term | Description |
|------|-------------|
| **npixmin** | Minimum number of pixels for a valid peak |
| **nsigma** | Significance threshold relative to background noise |
| **logrange** | Dynamic range for logarithmic contour levels |
| **logspacing** | Spacing between logarithmic contour levels |

### Iteration Parameters

| Term | Description |
|------|-------------|
| **iter_criterion** | Fractional change in λ required for convergence |
| **iter_nmax** | Maximum number of iterations |
| **initial_guess** | Starting estimate for λ in parsecs |

## Algorithms

### CLFIND (Clumpfind)

Hierarchical clump-finding algorithm (Williams, de Geus & Blitz 1994) used for peak identification:

1. Process contour levels from highest to lowest
2. At each level, identify connected regions (8-connectivity)
3. New peaks at local maxima, existing peaks grow downward
4. When clumps merge, assign pixels to nearest peak
5. Filter clumps by minimum pixel count (npixmin)

### KL14 Model

The analytical model from Kruijssen & Longmore (2014) that relates the flux ratio variance to the underlying timescales. Key functions:

| Function | Meaning |
|----------|---------|
| **f_zeta** | Dimensionless ratio function: (t_gas·t_star) / (t_total·t_over) |
| **f_rpeak** | Normalized peak radius function |
| **ζ (zeta)** | Dimensionless aperture size: l_ap / λ |

## Units

| Symbol | Unit | SI Value |
|--------|------|----------|
| pc | parsec | 3.086 × 10¹⁶ m |
| Myr | megayear | 3.156 × 10¹³ s |
| M☉ | solar mass | 1.989 × 10³⁰ kg |
| km/s | kilometers per second | 1000 m/s |

## Acronyms

| Acronym | Meaning |
|---------|---------|
| **KL14** | Kruijssen & Longmore (2014), the foundational paper |
| **GMC** | Giant Molecular Cloud |
| **SFR** | Star Formation Rate |
| **WCS** | World Coordinate System (astrometric information in FITS headers) |
| **DS9** | SAOImage DS9, astronomical image viewer |
| **TUI** | Text User Interface (interactive terminal mode) |
| **FITS** | Flexible Image Transport System (standard astronomical file format) |
| **FFT** | Fast Fourier Transform |
| **MC** | Monte Carlo (random sampling for error estimation) |

## References

- Kruijssen, J.M.D. & Longmore, S.N. (2014), MNRAS 439, 3239 - Original KL14 method paper
- Kruijssen et al. (2018), MNRAS 479, 1866 - Code paper and applications
- Williams, J.P., de Geus, E.J. & Blitz, L. (1994), ApJ 428, 693 - CLFIND algorithm

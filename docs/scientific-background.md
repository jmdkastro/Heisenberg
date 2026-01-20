# Scientific Background

Understanding the KL14 uncertainty principle for star formation.

## The Problem

Star formation is a multi-phase process:
1. Gas assembles into molecular clouds
2. Dense regions collapse to form stars
3. Stellar feedback disperses the remaining gas

**Key question:** How long does each phase last?

Traditional methods struggle because:
- We observe snapshots, not time evolution
- Individual clouds are hard to track
- Timescales vary from region to region

## The KL14 Uncertainty Principle

Kruijssen & Longmore (2014) developed a statistical method that exploits the **spatial decorrelation** between gas and young stellar tracers to derive timescales.

### Core Insight

On **large scales** (many star-forming regions):
- Gas and stars are spatially correlated
- Averaging over many regions gives the galactic mean

On **small scales** (individual regions):
- Gas and stars are decorrelated
- Each region is in a different evolutionary phase
- Some show only gas, some only stars, some both

The transition scale between these regimes encodes information about the underlying timescales.

### The Tuning Fork Diagram

When plotting the gas-to-stellar flux ratio as a function of aperture size:

```
Flux ratio
    ^
    |     *
    |    * *
    |   *   *     ← Small scales: large scatter
    |  *     *       (individual regions resolved)
    | *       *
    |*---------*--→ Large scales: converges to mean
    +------------→ Aperture size
```

The characteristic "tuning fork" shape arises because:
- Small apertures catch regions in different phases
- Large apertures average over many regions
- The two branches (gas-rich and star-rich) merge at scales > λ

## Physical Model

### Timeline Definitions

| Parameter | Symbol | Definition |
|-----------|--------|------------|
| Gas lifetime | t_gas | Time from cloud assembly to complete dispersal |
| Stellar timescale | t_star | Visibility timescale of stellar tracer |
| Overlap time | t_over | Duration when both tracers visible |
| Total timeline | t_total | t_gas + t_star - t_over |

### The Evolutionary Sequence

```
Time →

Phase 1: Gas only (t_gas - t_over)
|████████████████████████████|
         ↓
Phase 2: Overlap (t_over)
|████████████████████████████|
|▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓|  ← Both tracers visible
         ↓
Phase 3: Stars only (t_star - t_over)
|▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓|

█ = Gas tracer
▓ = Stellar tracer
```

### Region Separation Length (λ)

λ is the characteristic spacing between independent star-forming regions:
- Regions separated by > λ evolve independently
- Regions closer than λ may be correlated

λ relates to the physics of the ISM:
- Toomre length in galactic discs
- Jeans length in dense gas
- Feedback-regulated spacing

## The Fitting Procedure

### 1. Peak Identification

Identify emission peaks in both maps using the CLFIND algorithm:
- Contour-based hierarchical decomposition
- Each peak represents a star-forming region
- Characterize peak positions and fluxes

### 2. Aperture Photometry

For each peak, measure enclosed flux at multiple aperture sizes:
- Logarithmically spaced from lapmin to lapmax
- Typically 7-15 aperture sizes

### 3. Monte Carlo Sampling

Randomly sample peaks and compute flux ratios:
- Sample gas peaks and stellar peaks
- Calculate flux ratio at each aperture
- Repeat many times (nmc ~ 1000)
- Build probability distributions

### 4. Model Fitting

Fit the KL14 model to the observed tuning fork:
- Three free parameters: t_gas, t_over, λ
- One input parameter: t_star (from tracer calibration)
- Minimize chi-squared on a parameter grid
- Refine with adaptive grid search

### 5. Derived Quantities

From the fitted parameters, derive:
- Total timeline: t_total = t_gas + t_star - t_over
- Phase fractions
- Star formation efficiency (with conversion factors)
- Feedback velocity
- Mass loading factor

## Iterative Diffuse Filtering

Extended ("diffuse") emission can bias the analysis. The iterative procedure:

1. Run initial fit to get λ₀
2. Filter emission at scales > λ₀ using Fourier methods
3. Re-identify peaks in filtered maps
4. Refit to get λ₁
5. Repeat until |λₙ - λₙ₋₁|/λₙ < criterion

This converges because:
- Filtering removes large-scale emission
- Peaks become better defined
- λ estimate improves
- Filtering scale adjusts accordingly

## Assumptions and Limitations

### Key Assumptions

1. **Random sampling**: Regions observed are a random sample of the evolutionary sequence
2. **Constant timescales**: All regions follow the same timeline (population average)
3. **Independent regions**: Regions separated by > λ evolve independently
4. **Complete sampling**: Peak detection captures the full population

### Limitations

1. **Resolution**: Cannot probe scales smaller than the beam
2. **Statistics**: Requires sufficient number of peaks (~30+)
3. **Tracer timescale**: t_star must be known from other methods
4. **Stochasticity**: Single regions may deviate from population mean

## Applications

Heisenberg has been applied to:

- **Nearby galaxies**: M31, M33, NGC 300, Local Group dwarfs
- **PHANGS survey**: ~80 nearby star-forming galaxies
- **High-redshift**: Adapted methods for z ~ 1-3 galaxies
- **Simulations**: Validation against known input parameters

Key results include:
- GMC lifetimes of 10-30 Myr
- Short feedback timescales (1-5 Myr)
- Rapid, inefficient star formation
- Environmental dependence of timescales

## References

### Foundational Papers

- **Kruijssen & Longmore (2014)**, MNRAS 439, 3239
  - Original method paper
  - Derivation of the uncertainty principle

- **Kruijssen et al. (2018)**, MNRAS 479, 1866
  - Code paper
  - Detailed methodology
  - Application to M33

### Tracer Timescales

- **Haydon et al. (2020)**, MNRAS 498, 235
  - Calibration of stellar tracer timescales
  - Age dating with different tracers

### Applications

- **Chevance et al. (2020)**, MNRAS 493, 2872
  - PHANGS-ALMA survey results
  - 9 nearby galaxies

- **Kim et al. (2022)**, MNRAS 516, 3006
  - Cloud-scale star formation
  - Environmental variations

## Further Reading

For mathematical derivations and technical details, see:
- Kruijssen & Longmore (2014), Section 2
- Kruijssen et al. (2018), Section 3

For validation against simulations:
- Kruijssen et al. (2018), Section 5
- Hygate et al. (2019), MNRAS 488, 2800

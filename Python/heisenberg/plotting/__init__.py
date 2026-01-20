"""
Plotting module for Heisenberg.

Contains:
- Tuning fork diagnostic plots
- Image display functions (astronomical images)
- FITS file visualization
- Clump/peak visualization
- DS9 region file I/O
- Iteration history plots
"""

# Tuning fork plots
from heisenberg.plotting.tuningfork_plot import (
    plot_tuningfork,
    plot_tuningfork_residuals,
    plot_parameter_pdfs,
    plot_summary,
    save_figure,
)

# Image display (using astropy.visualization)
from heisenberg.plotting.image_display import (
    display_image,
    display_image_with_contours,
    get_percentile_scaling,
)

# FITS file display
from heisenberg.plotting.fits_display import (
    plot_fits,
    plot_fits_from_file,
    plot_fits_grid,
    plot_fits_grid_from_files,
    plot_star_gas_maps,
)

# Iteration history plots
from heisenberg.plotting.iteration_plot import (
    plot_iteration_history,
    plot_multiple_parameters,
    plot_convergence,
)

# DS9 region file I/O
from heisenberg.plotting.ds9_regions import (
    write_peak_regions,
    write_circle_regions,
    read_peak_regions,
    box_vertices,
    write_box_regions,
    read_box_regions,
)

# Clump visualization
from heisenberg.plotting.clump_display import (
    plot_clumps,
    plot_clumps_with_contours,
    plot_assignment_map,
    plot_clump_gallery,
    plot_peak_comparison,
)

__all__ = [
    # Tuning fork
    "plot_tuningfork",
    "plot_tuningfork_residuals",
    "plot_parameter_pdfs",
    "plot_summary",
    "save_figure",
    # Image display
    "display_image",
    "display_image_with_contours",
    "get_percentile_scaling",
    # FITS display
    "plot_fits",
    "plot_fits_from_file",
    "plot_fits_grid",
    "plot_fits_grid_from_files",
    "plot_star_gas_maps",
    # Iteration plots
    "plot_iteration_history",
    "plot_multiple_parameters",
    "plot_convergence",
    # DS9 regions
    "write_peak_regions",
    "write_circle_regions",
    "read_peak_regions",
    "box_vertices",
    "write_box_regions",
    "read_box_regions",
    # Clump display
    "plot_clumps",
    "plot_clumps_with_contours",
    "plot_assignment_map",
    "plot_clump_gallery",
    "plot_peak_comparison",
]

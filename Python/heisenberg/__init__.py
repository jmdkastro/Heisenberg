"""
Heisenberg: Uncertainty Principle for Star Formation

A Python implementation of the Heisenberg methodology for measuring GMC lifetime,
gas clearance timescale, and star formation efficiency from galaxy tracer maps.

Based on the methodology developed by Kruijssen & Longmore (2014) and
Kruijssen et al. (2018), with extensions by Hygate et al. (2019) and
Haydon et al. (2020a, b).

Example usage:
    >>> from heisenberg import HeisenbergConfig, run_analysis
    >>> config = HeisenbergConfig.from_file("input_file")
    >>> results = run_analysis(config)
"""

__version__ = "0.1.0"
__author__ = "Heisenberg Team"

from heisenberg.config import HeisenbergConfig, load_config
from heisenberg.core.constants import ASTRO, PHYS, NUMBERS

__all__ = [
    "__version__",
    "HeisenbergConfig",
    "load_config",
    "ASTRO",
    "PHYS",
    "NUMBERS",
]

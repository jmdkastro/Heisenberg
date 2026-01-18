"""
Configuration module for Heisenberg.

Provides the HeisenbergConfig dataclass and utilities for loading
configuration from IDL-format input files.

Example usage:
    >>> from heisenberg.config import HeisenbergConfig, load_config
    >>> config = load_config("path/to/input_file")
    >>> print(config.galaxy)
    >>> print(config.distance)
"""

from heisenberg.config.parameters import HeisenbergConfig
from heisenberg.config.reader import load_config, parse_input_file

__all__ = ["HeisenbergConfig", "load_config", "parse_input_file"]

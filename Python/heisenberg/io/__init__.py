"""
I/O module for Heisenberg.

Contains FITS file reading/writing and output file generation.
"""

from heisenberg.io.output import (
    OutputConfig,
    write_results,
    read_results,
    write_json,
    read_json,
    write_yaml,
    read_yaml,
    write_summary,
    result_to_dict,
)

__all__ = [
    "OutputConfig",
    "write_results",
    "read_results",
    "write_json",
    "read_json",
    "write_yaml",
    "read_yaml",
    "write_summary",
    "result_to_dict",
]

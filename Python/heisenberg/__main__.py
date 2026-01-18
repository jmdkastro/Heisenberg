"""
Entry point for running Heisenberg as a module.

Usage:
    python -m heisenberg --help
    python -m heisenberg run input_file
    python -m heisenberg run input_file --no-diffuse-filtering
"""

from heisenberg.cli import main

if __name__ == "__main__":
    main()

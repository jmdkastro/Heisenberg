"""
Physical and astronomical constants for the Heisenberg analysis.

These constants are extracted from the IDL COMMON blocks in heisenberg_nodf.pro.
All values are in SI units unless otherwise noted.

Constants are organized into three frozen dataclasses:
- Numbers: Numerical limits and dynamic range
- PhysConst: Fundamental physical constants
- AstrConst: Astronomical constants and unit conversions

Example usage:
    >>> from heisenberg.core.constants import ASTRO, PHYS, NUMBERS
    >>> distance_m = 10.0 * ASTRO.kpc  # 10 kpc in meters
    >>> print(f"Distance: {distance_m:.2e} m")
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Numbers:
    """
    Numerical limits and dynamic range for computations.

    Attributes:
        tiny: Very small number to avoid division by zero (1e-99)
        huge: Very large number for upper bounds (9e99)
        dynrange: Dynamic range for numerical operations (1e10)
    """
    tiny: float = 1.0e-99
    huge: float = 9.0e99
    dynrange: float = 1.0e10


@dataclass(frozen=True)
class PhysConst:
    """
    Fundamental physical constants in SI units.

    Attributes:
        kboltz: Boltzmann constant [J K^-1]
        ggravity: Gravitational constant [N m^2 kg^-2]
        hplanck: Planck constant [J s]
        clight: Speed of light in vacuum [m s^-1]
        mhydrogen: Hydrogen atom mass [kg]
        sigmaboltz: Stefan-Boltzmann constant [W m^-2 K^-4]
    """
    kboltz: float = 1.3807e-23       # Boltzmann constant [J K^-1]
    ggravity: float = 6.67259e-11    # Gravitational constant [N m^2 kg^-2]
    hplanck: float = 6.62607e-34     # Planck constant [J s]
    clight: float = 2.99792458e8     # Speed of light [m s^-1]
    mhydrogen: float = 1.6726e-27    # Hydrogen mass [kg]
    sigmaboltz: float = 5.67e-8      # Stefan-Boltzmann constant [W m^-2 K^-4]


@dataclass(frozen=True)
class AstrConst:
    """
    Astronomical constants and unit conversions in SI units.

    Attributes:
        kpc: Kiloparsec in meters [m]
        pc: Parsec in meters [m]
        msun: Solar mass [kg]
        lsun: Solar luminosity [W]
        yr: Sidereal year [s]
        myr: Million years [s]
        kms: km/s in m/s [m s^-1]
    """
    kpc: float = 3.08572e19          # Kiloparsec [m]
    pc: float = 3.08572e16           # Parsec [m] (kpc / 1000)
    msun: float = 1.989e30           # Solar mass [kg]
    lsun: float = 3.862e26           # Solar luminosity [W]
    yr: float = 86400.0 * 365.25     # Sidereal year [s]
    myr: float = 1.0e6 * 86400.0 * 365.25  # Million years [s]
    kms: float = 1.0e3               # km/s to m/s [m s^-1]

    @property
    def gyr(self) -> float:
        """Billion years in seconds."""
        return self.myr * 1.0e3


# Module-level singleton instances for easy import
NUMBERS = Numbers()
PHYS = PhysConst()
ASTRO = AstrConst()


# Convenience aliases matching common usage
tiny = NUMBERS.tiny
huge = NUMBERS.huge
dynrange = NUMBERS.dynrange

kboltz = PHYS.kboltz
ggravity = PHYS.ggravity
hplanck = PHYS.hplanck
clight = PHYS.clight
mhydrogen = PHYS.mhydrogen
sigmaboltz = PHYS.sigmaboltz

kpc = ASTRO.kpc
pc = ASTRO.pc
msun = ASTRO.msun
lsun = ASTRO.lsun
yr = ASTRO.yr
myr = ASTRO.myr
kms = ASTRO.kms

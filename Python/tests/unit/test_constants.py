"""
Tests for physical and astronomical constants.
"""

import pytest
from heisenberg.core.constants import NUMBERS, PHYS, ASTRO


class TestNumbers:
    """Tests for numerical constants."""

    def test_tiny_is_very_small(self):
        """Tiny should be very small to avoid division by zero."""
        assert NUMBERS.tiny == 1.0e-99
        assert NUMBERS.tiny > 0

    def test_huge_is_very_large(self):
        """Huge should be very large for upper bounds."""
        assert NUMBERS.huge == 9.0e99

    def test_dynrange(self):
        """Dynamic range should be 1e10."""
        assert NUMBERS.dynrange == 1.0e10


class TestPhysConst:
    """Tests for physical constants."""

    def test_boltzmann(self):
        """Boltzmann constant should be ~1.38e-23 J/K."""
        assert PHYS.kboltz == pytest.approx(1.3807e-23, rel=1e-4)

    def test_gravitational(self):
        """Gravitational constant should be ~6.67e-11 N m^2 kg^-2."""
        assert PHYS.ggravity == pytest.approx(6.67259e-11, rel=1e-4)

    def test_planck(self):
        """Planck constant should be ~6.63e-34 J s."""
        assert PHYS.hplanck == pytest.approx(6.62607e-34, rel=1e-4)

    def test_speed_of_light(self):
        """Speed of light should be ~3e8 m/s."""
        assert PHYS.clight == pytest.approx(2.99792458e8, rel=1e-9)

    def test_hydrogen_mass(self):
        """Hydrogen mass should be ~1.67e-27 kg."""
        assert PHYS.mhydrogen == pytest.approx(1.6726e-27, rel=1e-4)


class TestAstrConst:
    """Tests for astronomical constants."""

    def test_parsec(self):
        """Parsec should be ~3.086e16 m."""
        assert ASTRO.pc == pytest.approx(3.08572e16, rel=1e-4)

    def test_kiloparsec(self):
        """Kiloparsec should be 1000 parsecs."""
        assert ASTRO.kpc == pytest.approx(1000 * ASTRO.pc, rel=1e-10)

    def test_solar_mass(self):
        """Solar mass should be ~1.99e30 kg."""
        assert ASTRO.msun == pytest.approx(1.989e30, rel=1e-4)

    def test_solar_luminosity(self):
        """Solar luminosity should be ~3.86e26 W."""
        assert ASTRO.lsun == pytest.approx(3.862e26, rel=1e-4)

    def test_year_seconds(self):
        """Year should be ~3.16e7 seconds."""
        assert ASTRO.yr == pytest.approx(3.15576e7, rel=1e-4)

    def test_megayear(self):
        """Megayear should be 1e6 years."""
        assert ASTRO.myr == pytest.approx(1e6 * ASTRO.yr, rel=1e-10)

    def test_gigayear_property(self):
        """Gigayear should be 1e3 megayears."""
        assert ASTRO.gyr == pytest.approx(1e3 * ASTRO.myr, rel=1e-10)

    def test_kms(self):
        """km/s conversion should be 1000."""
        assert ASTRO.kms == 1000.0


class TestConstantImmutability:
    """Tests that constants are immutable."""

    def test_numbers_frozen(self):
        """Numbers dataclass should be frozen."""
        with pytest.raises(Exception):  # FrozenInstanceError
            NUMBERS.tiny = 1.0

    def test_phys_frozen(self):
        """PhysConst dataclass should be frozen."""
        with pytest.raises(Exception):
            PHYS.clight = 1.0

    def test_astro_frozen(self):
        """AstrConst dataclass should be frozen."""
        with pytest.raises(Exception):
            ASTRO.pc = 1.0

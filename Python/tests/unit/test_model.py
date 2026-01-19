"""Tests for heisenberg.core.model functions."""

import numpy as np
import pytest

from heisenberg.core.model import (
    f_zeta,
    f_rpeak,
    f_fluxratiostar,
    f_fluxratiogas,
    gaussian,
)
from heisenberg.core.constants import NUMBERS


class TestFZeta:
    """Tests for f_zeta function."""

    def test_point_profile_returns_tiny(self):
        """Peak profile 0 (point) should return tiny."""
        result = f_zeta(surfcontrast=10.0, ttot=10.0, tref=5.0, peak_prof=0)
        assert result == NUMBERS.tiny

    def test_disc_profile(self):
        """Peak profile 1 (disc) calculation."""
        # With surfcontrast=2, ttot=10, tref=5:
        # zeta = sqrt((10/5) / (2*(1+10/5) - 1)) = sqrt(2 / (2*3 - 1)) = sqrt(2/5)
        result = f_zeta(surfcontrast=2.0, ttot=10.0, tref=5.0, peak_prof=1)
        expected = np.sqrt(2.0 / 5.0)
        assert np.isclose(result, expected, rtol=1e-10)

    def test_gaussian_profile(self):
        """Peak profile 2 (Gaussian) calculation."""
        # With surfcontrast=2, ttot=10, tref=5:
        # zeta = sqrt((10/5) / (2*2*(1+10/5) - 2)) = sqrt(2 / (4*3 - 2)) = sqrt(2/10)
        result = f_zeta(surfcontrast=2.0, ttot=10.0, tref=5.0, peak_prof=2)
        expected = np.sqrt(2.0 / 10.0)
        assert np.isclose(result, expected, rtol=1e-10)

    def test_surfcontrast_clipped_to_one(self):
        """Surface contrast below 1 should be clipped to 1."""
        result = f_zeta(surfcontrast=0.5, ttot=10.0, tref=5.0, peak_prof=1)
        result_at_one = f_zeta(surfcontrast=1.0, ttot=10.0, tref=5.0, peak_prof=1)
        assert result == result_at_one

    def test_invalid_peak_prof_raises(self):
        """Invalid peak_prof should raise ValueError."""
        with pytest.raises(ValueError):
            f_zeta(surfcontrast=2.0, ttot=10.0, tref=5.0, peak_prof=3)


class TestFRpeak:
    """Tests for f_rpeak function."""

    def test_basic_calculation(self):
        """Basic rpeak = 0.5 * zeta * lambda."""
        result = f_rpeak(zeta=0.5, lambda_=100.0)
        assert result == 25.0

    def test_zero_zeta(self):
        """Zero zeta gives zero rpeak."""
        result = f_rpeak(zeta=0.0, lambda_=100.0)
        assert result == 0.0


class TestFluxRatioStar:
    """Tests for f_fluxratiostar function."""

    def test_output_shape(self):
        """Output should match aperture array shape."""
        apertures = np.array([50, 100, 200, 400])
        result = f_fluxratiostar(
            tgas=10.0,
            tstar=5.0,
            tover=1.0,
            laps=apertures,
            lambda_=100.0,
            beta_gas=1.0,
            surfcontrasts=5.0,
            surfcontrastg=5.0,
            peak_prof=2
        )
        assert result.shape == apertures.shape

    def test_flux_ratio_positive(self):
        """Flux ratios should be positive."""
        apertures = np.array([50, 100, 200, 400])
        result = f_fluxratiostar(
            tgas=10.0,
            tstar=5.0,
            tover=1.0,
            laps=apertures,
            lambda_=100.0,
            beta_gas=1.0,
            surfcontrasts=5.0,
            surfcontrastg=5.0,
            peak_prof=2
        )
        assert np.all(result > 0)

    def test_galactic_average_limit(self):
        """At very large apertures, flux ratio should approach 1."""
        apertures = np.array([10000, 50000, 100000])
        result = f_fluxratiostar(
            tgas=10.0,
            tstar=5.0,
            tover=1.0,
            laps=apertures,
            lambda_=100.0,
            beta_gas=1.0,
            surfcontrasts=5.0,
            surfcontrastg=5.0,
            peak_prof=2
        )
        # At large apertures, should approach 1 (galactic average)
        assert np.all(np.abs(result - 1.0) < 0.1)


class TestFluxRatioGas:
    """Tests for f_fluxratiogas function."""

    def test_output_shape(self):
        """Output should match aperture array shape."""
        apertures = np.array([50, 100, 200, 400])
        result = f_fluxratiogas(
            tgas=10.0,
            tstar=5.0,
            tover=1.0,
            lapg=apertures,
            lambda_=100.0,
            beta_star=1.0,
            surfcontrasts=5.0,
            surfcontrastg=5.0,
            peak_prof=2
        )
        assert result.shape == apertures.shape

    def test_flux_ratio_positive(self):
        """Flux ratios should be positive."""
        apertures = np.array([50, 100, 200, 400])
        result = f_fluxratiogas(
            tgas=10.0,
            tstar=5.0,
            tover=1.0,
            lapg=apertures,
            lambda_=100.0,
            beta_star=1.0,
            surfcontrasts=5.0,
            surfcontrastg=5.0,
            peak_prof=2
        )
        assert np.all(result > 0)

    def test_gas_exceeds_star_at_small_apertures(self):
        """Gas flux ratio should exceed stellar at small apertures (tuning fork)."""
        apertures = np.array([50, 100, 200])
        result_gas = f_fluxratiogas(
            tgas=10.0,
            tstar=5.0,
            tover=1.0,
            lapg=apertures,
            lambda_=100.0,
            beta_star=1.0,
            surfcontrasts=5.0,
            surfcontrastg=5.0,
            peak_prof=2
        )
        result_star = f_fluxratiostar(
            tgas=10.0,
            tstar=5.0,
            tover=1.0,
            laps=apertures,
            lambda_=100.0,
            beta_gas=1.0,
            surfcontrasts=5.0,
            surfcontrastg=5.0,
            peak_prof=2
        )
        # Gas should show higher depletion time (more gas per star)
        # Star should show lower depletion time (more star per gas)
        # So gas flux ratio > 1 and star flux ratio < 1 at small apertures
        assert result_gas[0] > result_star[0]


class TestGaussian:
    """Tests for gaussian function."""

    def test_peak_value(self):
        """Gaussian peak should equal amplitude."""
        params = np.array([5.0, 0.0, 1.0])  # amplitude, center, sigma
        result = gaussian(np.array([0.0]), params)
        assert np.isclose(result[0], 5.0)

    def test_symmetric(self):
        """Gaussian should be symmetric around center."""
        params = np.array([1.0, 2.0, 1.0])
        x = np.array([1.0, 3.0])  # Equidistant from center
        result = gaussian(x, params)
        assert np.isclose(result[0], result[1])

    def test_fwhm(self):
        """Check FWHM relationship: FWHM = 2*sqrt(2*ln(2))*sigma."""
        sigma = 1.0
        params = np.array([1.0, 0.0, sigma])
        fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma
        half_max_x = fwhm / 2.0
        result = gaussian(np.array([half_max_x]), params)
        assert np.isclose(result[0], 0.5, rtol=1e-10)

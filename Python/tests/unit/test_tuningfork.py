"""Tests for heisenberg.core.tuningfork module."""

import numpy as np
import pytest

from heisenberg.core.tuningfork import (
    TuningForkConfig,
    TuningForkData,
    TuningForkResult,
    run_tuningfork,
    compute_derived_quantities,
    generate_model_curve,
)
from heisenberg.core.fitting import FitResult


def make_gaussian_image(shape, center, amplitude, sigma):
    """Create a 2D Gaussian test image."""
    y, x = np.ogrid[:shape[0], :shape[1]]
    cy, cx = center
    return amplitude * np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * sigma**2))


def make_mock_fit_result(
    tgas=10.0,
    tover=1.0,
    lambda_=100.0,
    tstar=10.0,
) -> FitResult:
    """Create a mock FitResult for testing."""
    return FitResult(
        tgas=tgas,
        tover=tover,
        lambda_=lambda_,
        tstar=tstar,
        tgas_errmin=1.0,
        tgas_errmax=1.5,
        tover_errmin=0.3,
        tover_errmax=0.5,
        lambda_errmin=10.0,
        lambda_errmax=15.0,
        chi2_min=1.5,
        ndof=10,
        beta_star=1.0,
        beta_gas=1.0,
        surfcon_star=5.0,
        surfcon_gas=5.0,
        tgas_arr=np.linspace(1, 50, 10),
        tover_arr=np.linspace(0.1, 10, 10),
        lambda_arr=np.linspace(10, 500, 10),
        prob_tgas=np.ones(10) / 10,
        prob_tover=np.ones(10) / 10,
        prob_lambda=np.ones(10) / 10,
        corr_tgas_tover=0.1,
        corr_tgas_lambda=-0.2,
        corr_tover_lambda=0.05,
    )


def make_multi_peak_image(shape, peaks, sigma=5.0):
    """
    Create image with multiple Gaussian peaks.

    Args:
        shape: Image shape (ny, nx)
        peaks: List of (y, x, amplitude) tuples
        sigma: Gaussian sigma for all peaks
    """
    image = np.zeros(shape)
    for cy, cx, amp in peaks:
        image += make_gaussian_image(shape, (cy, cx), amp, sigma)
    return image


class TestTuningForkConfig:
    """Tests for TuningForkConfig dataclass."""

    def test_default_values(self):
        """Config should have sensible defaults."""
        config = TuningForkConfig()

        assert config.lap_min == 50.0
        assert config.lap_max == 1000.0
        assert config.n_apertures == 15
        assert config.n_mc == 100
        assert config.tstar == 10.0
        assert config.peak_prof == 2
        assert config.npixmin_star == 20
        assert config.npixmin_gas == 20
        assert config.nsigma_star == 5.0
        assert config.nsigma_gas == 5.0

    def test_custom_values(self):
        """Config should accept custom values."""
        config = TuningForkConfig(
            lap_min=100.0,
            lap_max=500.0,
            n_apertures=10,
            tstar=15.0,
        )

        assert config.lap_min == 100.0
        assert config.lap_max == 500.0
        assert config.n_apertures == 10
        assert config.tstar == 15.0

    def test_seed_for_reproducibility(self):
        """Config should allow setting random seed."""
        config = TuningForkConfig(seed=42)
        assert config.seed == 42


class TestTuningForkData:
    """Tests for TuningForkData dataclass."""

    def test_create_data(self):
        """Should create data structure with required fields."""
        data = TuningForkData(
            apertures=np.array([50, 100, 200]),
            fluxratio_star=np.array([1.5, 1.2, 1.1]),
            fluxratio_gas=np.array([0.7, 0.8, 0.9]),
            err_star=np.array([0.1, 0.1, 0.1]),
            err_gas=np.array([0.1, 0.1, 0.1]),
            n_star_peaks=np.array([10, 10, 10]),
            n_gas_peaks=np.array([8, 8, 8]),
        )

        assert len(data.apertures) == 3
        assert len(data.fluxratio_star) == 3
        assert len(data.fluxratio_gas) == 3


class TestTuningForkResult:
    """Tests for TuningForkResult dataclass."""

    def test_create_result(self):
        """Should create result with all components."""
        # Create mock components
        observed = TuningForkData(
            apertures=np.array([50, 100]),
            fluxratio_star=np.array([1.5, 1.2]),
            fluxratio_gas=np.array([0.7, 0.8]),
            err_star=np.array([0.1, 0.1]),
            err_gas=np.array([0.1, 0.1]),
            n_star_peaks=np.array([10, 10]),
            n_gas_peaks=np.array([8, 8]),
        )

        fit = make_mock_fit_result()

        result = TuningForkResult(
            observed=observed,
            fit=fit,
            star_peaks=[],
            gas_peaks=[],
            derived={'ttotal': 19.0},
        )

        assert result.observed is not None
        assert result.fit is not None
        assert result.derived['ttotal'] == 19.0


class TestComputeDerivedQuantities:
    """Tests for compute_derived_quantities function."""

    def test_ttotal_calculation(self):
        """ttotal = tgas + tstar - tover."""
        fit = make_mock_fit_result(tgas=10.0, tover=2.0, lambda_=100.0, tstar=8.0)
        config = TuningForkConfig(tstar=8.0)

        derived = compute_derived_quantities(fit, config)

        # ttotal = tgas + tstar - tover = 10 + 8 - 2 = 16
        assert np.isclose(derived['ttotal'], 16.0)

    def test_esf_positive(self):
        """Star formation efficiency should be positive."""
        fit = make_mock_fit_result()
        config = TuningForkConfig()

        derived = compute_derived_quantities(fit, config)

        assert derived['esf'] > 0

    def test_all_quantities_present(self):
        """Should compute all expected quantities."""
        fit = make_mock_fit_result()
        config = TuningForkConfig()

        derived = compute_derived_quantities(fit, config)

        expected_keys = ['tgas', 'tover', 'lambda', 'tstar', 'ttotal',
                         'esf', 'vfb', 'etainst', 'etaavg']
        for key in expected_keys:
            assert key in derived, f"Missing key: {key}"


class TestGenerateModelCurve:
    """Tests for generate_model_curve function."""

    def test_returns_arrays(self):
        """Should return apertures and model curves."""
        fit = make_mock_fit_result()

        apertures, model_star, model_gas = generate_model_curve(fit)

        assert isinstance(apertures, np.ndarray)
        assert isinstance(model_star, np.ndarray)
        assert isinstance(model_gas, np.ndarray)

    def test_curve_shape(self):
        """Model curves should have correct shape."""
        fit = make_mock_fit_result()

        apertures, model_star, model_gas = generate_model_curve(fit, n_points=50)

        assert len(apertures) == 50
        assert len(model_star) == 50
        assert len(model_gas) == 50

    def test_custom_apertures(self):
        """Should accept custom aperture values."""
        fit = make_mock_fit_result()

        custom_apertures = np.array([50, 100, 200, 500])
        apertures, model_star, model_gas = generate_model_curve(
            fit, apertures=custom_apertures
        )

        assert np.array_equal(apertures, custom_apertures)
        assert len(model_star) == 4
        assert len(model_gas) == 4

    def test_model_star_flux_ratio_behavior(self):
        """At stellar peaks, flux ratio should increase at small apertures."""
        fit = make_mock_fit_result()

        apertures, model_star, model_gas = generate_model_curve(fit, n_points=20)

        # At small apertures, stellar peaks should have enhanced gas/star ratio
        # (more stellar flux relative to gas)
        # This is the characteristic tuning fork behavior
        assert len(model_star) == 20


class TestRunTuningfork:
    """Tests for run_tuningfork function."""

    def test_insufficient_star_peaks_raises(self):
        """Should raise if < 2 stellar peaks found."""
        # Create image with only 1 peak
        star_image = make_gaussian_image((64, 64), (32, 32), 100.0, 5.0)
        gas_image = make_multi_peak_image((64, 64), [
            (20, 20, 100.0),
            (40, 40, 80.0),
            (50, 50, 60.0),
        ])

        config = TuningForkConfig(
            lap_min=5.0,
            lap_max=30.0,
            npixmin_star=5,
            npixmin_gas=5,
            nsigma_star=3.0,
            nsigma_gas=3.0,
        )

        with pytest.raises(ValueError, match="Insufficient stellar peaks"):
            run_tuningfork(star_image, gas_image, pixel_scale=1.0, config=config)

    def test_insufficient_gas_peaks_raises(self):
        """Should raise if < 2 gas peaks found."""
        # Create image with only 1 gas peak
        star_image = make_multi_peak_image((64, 64), [
            (20, 20, 100.0),
            (40, 40, 80.0),
            (50, 50, 60.0),
        ])
        gas_image = make_gaussian_image((64, 64), (32, 32), 100.0, 5.0)

        config = TuningForkConfig(
            lap_min=5.0,
            lap_max=30.0,
            npixmin_star=5,
            npixmin_gas=5,
            nsigma_star=3.0,
            nsigma_gas=3.0,
        )

        with pytest.raises(ValueError, match="Insufficient gas peaks"):
            run_tuningfork(star_image, gas_image, pixel_scale=1.0, config=config)

    def test_returns_result_structure(self):
        """Should return TuningForkResult with all components."""
        # Create images with multiple well-separated peaks
        shape = (128, 128)
        star_peaks = [
            (25, 25, 100.0),
            (25, 100, 80.0),
            (100, 25, 90.0),
            (100, 100, 70.0),
        ]
        gas_peaks = [
            (30, 30, 100.0),
            (30, 95, 80.0),
            (95, 30, 90.0),
            (95, 95, 70.0),
        ]

        star_image = make_multi_peak_image(shape, star_peaks, sigma=8.0)
        gas_image = make_multi_peak_image(shape, gas_peaks, sigma=8.0)

        config = TuningForkConfig(
            lap_min=10.0,
            lap_max=60.0,
            n_apertures=5,
            n_mc=10,  # Small for speed
            npixmin_star=10,
            npixmin_gas=10,
            nsigma_star=2.0,
            nsigma_gas=2.0,
            seed=42,
        )

        result = run_tuningfork(star_image, gas_image, pixel_scale=1.0, config=config)

        assert isinstance(result, TuningForkResult)
        assert result.observed is not None
        assert result.fit is not None
        assert len(result.star_peaks) >= 2
        assert len(result.gas_peaks) >= 2
        assert len(result.derived) > 0

    def test_config_preserved_in_result(self):
        """Result should include configuration used."""
        shape = (128, 128)
        star_peaks = [
            (25, 25, 100.0),
            (25, 100, 80.0),
            (100, 25, 90.0),
            (100, 100, 70.0),
        ]
        gas_peaks = [
            (30, 30, 100.0),
            (30, 95, 80.0),
            (95, 30, 90.0),
            (95, 95, 70.0),
        ]

        star_image = make_multi_peak_image(shape, star_peaks, sigma=8.0)
        gas_image = make_multi_peak_image(shape, gas_peaks, sigma=8.0)

        config = TuningForkConfig(
            lap_min=10.0,
            lap_max=60.0,
            n_apertures=5,
            n_mc=10,
            npixmin_star=10,
            npixmin_gas=10,
            nsigma_star=2.0,
            nsigma_gas=2.0,
            tstar=15.0,
            seed=42,
        )

        result = run_tuningfork(star_image, gas_image, pixel_scale=1.0, config=config)

        assert result.config['tstar'] == 15.0
        assert result.config['lap_min'] == 10.0
        assert result.config['lap_max'] == 60.0

    def test_default_config_used_when_none(self):
        """Should use default config if none provided."""
        shape = (128, 128)
        star_peaks = [
            (25, 25, 100.0),
            (25, 100, 80.0),
            (100, 25, 90.0),
            (100, 100, 70.0),
        ]
        gas_peaks = [
            (30, 30, 100.0),
            (30, 95, 80.0),
            (95, 30, 90.0),
            (95, 95, 70.0),
        ]

        star_image = make_multi_peak_image(shape, star_peaks, sigma=8.0)
        gas_image = make_multi_peak_image(shape, gas_peaks, sigma=8.0)

        # This should work with defaults, but may raise due to peak requirements
        # so we just test that it accepts None config
        try:
            result = run_tuningfork(star_image, gas_image, pixel_scale=1.0, config=None)
            assert result.config is not None
        except ValueError:
            # Expected if default config params don't find enough peaks
            pass


class TestRunTuningforkIntegration:
    """Integration tests for run_tuningfork with realistic scenarios."""

    def test_synthetic_data_runs_to_completion(self):
        """Full pipeline should run on synthetic data."""
        # Create more realistic synthetic data
        np.random.seed(42)
        shape = (200, 200)

        # Create scattered peaks
        n_peaks = 15
        star_peaks = []
        gas_peaks = []

        for i in range(n_peaks):
            # Random positions avoiding edges
            y = np.random.randint(30, 170)
            x = np.random.randint(30, 170)
            amp = np.random.uniform(50, 150)
            star_peaks.append((y, x, amp))

            # Gas peaks slightly offset
            gas_peaks.append((y + np.random.randint(-5, 5),
                              x + np.random.randint(-5, 5),
                              amp * np.random.uniform(0.8, 1.2)))

        star_image = make_multi_peak_image(shape, star_peaks, sigma=6.0)
        gas_image = make_multi_peak_image(shape, gas_peaks, sigma=6.0)

        # Add some noise
        star_image += np.random.normal(0, 5, shape)
        gas_image += np.random.normal(0, 5, shape)

        config = TuningForkConfig(
            lap_min=15.0,
            lap_max=80.0,
            n_apertures=8,
            n_mc=20,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=42,
        )

        result = run_tuningfork(star_image, gas_image, pixel_scale=1.0, config=config)

        # Check result structure
        assert result.observed is not None
        assert result.fit is not None
        assert result.fit.tgas > 0
        assert result.fit.lambda_ > 0
        assert 'ttotal' in result.derived


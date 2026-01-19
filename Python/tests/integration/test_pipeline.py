"""
Integration tests for the full Heisenberg analysis pipeline.

These tests verify that the complete pipeline works end-to-end
using synthetic data with known parameters.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from heisenberg.core.tuningfork import (
    TuningForkConfig,
    TuningForkResult,
    run_tuningfork,
    generate_model_curve,
)
from heisenberg.io import write_results, read_results, write_summary


def make_gaussian_image(shape, center, amplitude, sigma):
    """Create a 2D Gaussian test image."""
    y, x = np.ogrid[:shape[0], :shape[1]]
    cy, cx = center
    return amplitude * np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * sigma**2))


def make_multi_peak_image(shape, peaks, sigma=5.0, background=0.0):
    """
    Create image with multiple Gaussian peaks.

    Args:
        shape: Image shape (ny, nx)
        peaks: List of (y, x, amplitude) tuples
        sigma: Gaussian sigma for all peaks
        background: Background level to add
    """
    image = np.ones(shape) * background
    for cy, cx, amp in peaks:
        image += make_gaussian_image(shape, (cy, cx), amp, sigma)
    return image


def create_synthetic_galaxy(
    n_star_peaks: int = 20,
    n_gas_peaks: int = 20,
    shape: tuple = (256, 256),
    sigma: float = 8.0,
    seed: int = 42,
    star_flux_range: tuple = (50.0, 150.0),
    gas_flux_range: tuple = (40.0, 120.0),
    offset_range: tuple = (-5, 5),
    noise_level: float = 5.0,
):
    """
    Create synthetic star and gas maps simulating a galaxy.

    Args:
        n_star_peaks: Number of stellar peaks
        n_gas_peaks: Number of gas peaks
        shape: Image shape
        sigma: Gaussian sigma for peaks
        seed: Random seed
        star_flux_range: (min, max) amplitude for star peaks
        gas_flux_range: (min, max) amplitude for gas peaks
        offset_range: (min, max) offset for gas peaks from star peaks
        noise_level: Standard deviation of Gaussian noise

    Returns:
        Tuple of (star_image, gas_image)
    """
    np.random.seed(seed)

    margin = 30
    star_peaks = []
    gas_peaks = []

    for i in range(n_star_peaks):
        y = np.random.randint(margin, shape[0] - margin)
        x = np.random.randint(margin, shape[1] - margin)
        amp = np.random.uniform(*star_flux_range)
        star_peaks.append((y, x, amp))

    for i in range(n_gas_peaks):
        # Gas peaks slightly offset from star positions
        # Some overlap (correlated), some independent
        if i < n_star_peaks and np.random.random() < 0.7:
            # Correlated with stellar peak
            base_y, base_x, _ = star_peaks[i]
            y = base_y + np.random.randint(*offset_range)
            x = base_x + np.random.randint(*offset_range)
        else:
            # Independent gas peak
            y = np.random.randint(margin, shape[0] - margin)
            x = np.random.randint(margin, shape[1] - margin)

        amp = np.random.uniform(*gas_flux_range)
        gas_peaks.append((y, x, amp))

    star_image = make_multi_peak_image(shape, star_peaks, sigma)
    gas_image = make_multi_peak_image(shape, gas_peaks, sigma)

    # Add noise
    if noise_level > 0:
        star_image += np.random.normal(0, noise_level, shape)
        gas_image += np.random.normal(0, noise_level, shape)

    return star_image, gas_image


@pytest.mark.slow
class TestFullPipeline:
    """Integration tests for the complete analysis pipeline."""

    def test_pipeline_runs_to_completion(self):
        """Pipeline should complete without errors on synthetic data."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=15,
            n_gas_peaks=15,
            shape=(200, 200),
            seed=42,
        )

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

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        assert isinstance(result, TuningForkResult)
        assert result.fit is not None
        assert result.observed is not None
        assert len(result.star_peaks) >= 2
        assert len(result.gas_peaks) >= 2

    def test_pipeline_produces_positive_timescales(self):
        """Fitted timescales should be positive."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=20,
            n_gas_peaks=20,
            seed=123,
        )

        config = TuningForkConfig(
            lap_min=20.0,
            lap_max=100.0,
            n_apertures=6,
            n_mc=15,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=123,
        )

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        assert result.fit.tgas > 0
        assert result.fit.tover > 0
        assert result.fit.lambda_ > 0

    def test_pipeline_produces_physical_lambda(self):
        """Lambda should be positive and finite."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=15,
            n_gas_peaks=15,
            seed=456,
        )

        config = TuningForkConfig(
            lap_min=20.0,
            lap_max=100.0,
            n_apertures=6,
            n_mc=15,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=456,
        )

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        assert result.fit.lambda_ > 0
        assert np.isfinite(result.fit.lambda_)

    def test_derived_quantities_computed(self):
        """Pipeline should compute derived quantities."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=15,
            n_gas_peaks=15,
            seed=789,
        )

        config = TuningForkConfig(
            lap_min=20.0,
            lap_max=100.0,
            n_apertures=6,
            n_mc=15,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=789,
        )

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        assert 'ttotal' in result.derived
        assert 'esf' in result.derived
        assert 'vfb' in result.derived

    def test_model_curve_generation(self):
        """Should be able to generate model curves from fit."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=15,
            n_gas_peaks=15,
            seed=101,
        )

        config = TuningForkConfig(
            lap_min=20.0,
            lap_max=100.0,
            n_apertures=6,
            n_mc=15,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=101,
        )

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        apertures, model_star, model_gas = generate_model_curve(
            result.fit, n_points=50
        )

        assert len(apertures) == 50
        assert len(model_star) == 50
        assert len(model_gas) == 50
        assert np.all(np.isfinite(model_star))
        assert np.all(np.isfinite(model_gas))


@pytest.mark.slow
class TestPipelineWithMask:
    """Tests for pipeline with masked regions."""

    def test_pipeline_respects_mask(self):
        """Pipeline should work with masked regions."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=20,
            n_gas_peaks=20,
            seed=202,
        )

        # Create mask - exclude corners
        mask = np.ones(star_image.shape)
        mask[:50, :50] = 0  # Mask upper-left corner
        mask[-50:, -50:] = 0  # Mask lower-right corner

        config = TuningForkConfig(
            lap_min=20.0,
            lap_max=100.0,
            n_apertures=6,
            n_mc=15,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=202,
        )

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
            mask=mask,
        )

        assert result.fit is not None


@pytest.mark.slow
class TestPipelineIO:
    """Tests for pipeline output I/O integration."""

    def test_result_save_load_round_trip(self):
        """Results should survive save/load cycle."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=15,
            n_gas_peaks=15,
            seed=303,
        )

        config = TuningForkConfig(
            lap_min=20.0,
            lap_max=100.0,
            n_apertures=6,
            n_mc=15,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=303,
        )

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
            path = Path(f.name)

        try:
            write_results(result, path)
            loaded = read_results(path)

            # Compare key parameters
            assert np.isclose(loaded.fit.tgas, result.fit.tgas, rtol=1e-5)
            assert np.isclose(loaded.fit.tover, result.fit.tover, rtol=1e-5)
            assert np.isclose(loaded.fit.lambda_, result.fit.lambda_, rtol=1e-5)
        finally:
            path.unlink()

    def test_summary_generation(self):
        """Should generate readable summary file."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=15,
            n_gas_peaks=15,
            seed=404,
        )

        config = TuningForkConfig(
            lap_min=20.0,
            lap_max=100.0,
            n_apertures=6,
            n_mc=15,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=404,
        )

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            path = Path(f.name)

        try:
            write_summary(result, path)

            with open(path, 'r') as f:
                content = f.read()

            assert 't_gas' in content
            assert 't_over' in content
            assert 'lambda' in content
        finally:
            path.unlink()


@pytest.mark.slow
class TestPipelineReproducibility:
    """Tests for pipeline reproducibility with seeds."""

    def test_same_seed_same_result(self):
        """Same seed should produce same results."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=15,
            n_gas_peaks=15,
            seed=500,
        )

        config = TuningForkConfig(
            lap_min=20.0,
            lap_max=100.0,
            n_apertures=6,
            n_mc=15,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=500,
        )

        result1 = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        result2 = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        assert result1.fit.tgas == result2.fit.tgas
        assert result1.fit.tover == result2.fit.tover
        assert result1.fit.lambda_ == result2.fit.lambda_


@pytest.mark.slow
class TestPipelineEdgeCases:
    """Tests for edge cases and error handling."""

    def test_insufficient_peaks_raises(self):
        """Should raise clear error with insufficient peaks."""
        # Single peak only
        star_image = make_gaussian_image((100, 100), (50, 50), 100.0, 5.0)
        gas_image = make_multi_peak_image((100, 100), [
            (30, 30, 100.0),
            (70, 70, 80.0),
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

    def test_small_aperture_range(self):
        """Should handle narrow aperture range."""
        star_image, gas_image = create_synthetic_galaxy(
            n_star_peaks=15,
            n_gas_peaks=15,
            seed=600,
        )

        config = TuningForkConfig(
            lap_min=30.0,
            lap_max=50.0,  # Narrow range
            n_apertures=4,
            n_mc=10,
            npixmin_star=15,
            npixmin_gas=15,
            nsigma_star=3.0,
            nsigma_gas=3.0,
            seed=600,
        )

        result = run_tuningfork(
            star_image,
            gas_image,
            pixel_scale=1.0,
            config=config,
        )

        assert result.fit is not None

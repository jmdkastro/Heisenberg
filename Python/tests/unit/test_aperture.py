"""Tests for heisenberg.imaging.aperture module."""

import numpy as np
import pytest

from heisenberg.imaging.aperture import (
    ApertureFluxResult,
    generate_aperture_sizes,
    create_smoothed_cube,
    measure_aperture_flux,
    compute_flux_ratios,
    select_non_overlapping_peaks,
    monte_carlo_flux_ratios,
)


def make_gaussian_image(shape, center, amplitude, sigma):
    """Create a 2D Gaussian test image."""
    y, x = np.ogrid[:shape[0], :shape[1]]
    cy, cx = center
    return amplitude * np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * sigma**2))


class TestGenerateApertureSizes:
    """Tests for generate_aperture_sizes function."""

    def test_correct_count(self):
        """Should generate requested number of apertures."""
        apertures = generate_aperture_sizes(10, 100, 5)
        assert len(apertures) == 5

    def test_geometric_series(self):
        """Should be a geometric series."""
        apertures = generate_aperture_sizes(10, 100, 5)
        # Ratio between consecutive elements should be constant
        ratios = apertures[1:] / apertures[:-1]
        assert np.allclose(ratios, ratios[0])

    def test_min_max_bounds(self):
        """Should start at min and end at max."""
        apertures = generate_aperture_sizes(10, 100, 5)
        assert np.isclose(apertures[0], 10)
        assert np.isclose(apertures[-1], 100)

    def test_single_aperture(self):
        """Single aperture should return min value."""
        apertures = generate_aperture_sizes(10, 100, 1)
        assert len(apertures) == 1
        assert apertures[0] == 10


class TestCreateSmoothedCube:
    """Tests for create_smoothed_cube function."""

    def test_output_shape(self):
        """Output should have shape (n_apertures, ny, nx)."""
        image = np.random.randn(64, 64)
        apertures = np.array([5.0, 10.0, 20.0])
        cube = create_smoothed_cube(image, apertures, pixel_scale=1.0)
        assert cube.shape == (3, 64, 64)

    def test_increasing_smoothness(self):
        """Larger apertures should give smoother (lower variance) images."""
        np.random.seed(42)
        image = np.random.randn(64, 64)
        apertures = np.array([5.0, 10.0, 20.0])
        cube = create_smoothed_cube(image, apertures, pixel_scale=1.0)

        variances = [np.var(cube[i]) for i in range(3)]
        # Variance should decrease with larger apertures
        assert variances[0] >= variances[1] >= variances[2]


class TestMeasureApertureFlux:
    """Tests for measure_aperture_flux function."""

    def test_returns_two_results(self):
        """Should return results for both star and gas peaks."""
        image = make_gaussian_image((64, 64), (32, 32), 100.0, 5.0)
        star_peaks = np.array([[32, 32]])
        gas_peaks = np.array([[32, 32]])
        apertures = np.array([10.0, 20.0])

        star_result, gas_result = measure_aperture_flux(
            image, image, star_peaks, gas_peaks,
            apertures, pixel_scale=1.0
        )

        assert isinstance(star_result, ApertureFluxResult)
        assert isinstance(gas_result, ApertureFluxResult)

    def test_flux_array_shape(self):
        """Flux arrays should have shape (n_apertures, n_peaks)."""
        image = make_gaussian_image((64, 64), (32, 32), 100.0, 5.0)
        star_peaks = np.array([[32, 32], [40, 40]])
        gas_peaks = np.array([[32, 32]])
        apertures = np.array([10.0, 20.0, 30.0])

        star_result, gas_result = measure_aperture_flux(
            image, image, star_peaks, gas_peaks,
            apertures, pixel_scale=1.0
        )

        assert star_result.star_flux.shape == (3, 2)
        assert gas_result.gas_flux.shape == (3, 1)

    def test_flux_at_peak(self):
        """Flux at peak position should be positive."""
        image = make_gaussian_image((64, 64), (32, 32), 100.0, 5.0)
        peaks = np.array([[32, 32]])
        apertures = np.array([10.0])

        result, _ = measure_aperture_flux(
            image, image, peaks, peaks,
            apertures, pixel_scale=1.0
        )

        assert result.star_flux[0, 0] > 0


class TestComputeFluxRatios:
    """Tests for compute_flux_ratios function."""

    def test_returns_correct_shape(self):
        """Should return one ratio per aperture."""
        # Create mock flux result
        result = ApertureFluxResult(
            star_flux=np.array([[10.0, 20.0], [15.0, 25.0]]),
            gas_flux=np.array([[20.0, 40.0], [30.0, 50.0]]),
            aperture_area_frac=np.ones((2, 2)),
            apertures=np.array([10.0, 20.0]),
            peak_positions=np.array([[30, 30], [40, 40]]),
        )

        ratios, valid = compute_flux_ratios(result)

        assert len(ratios) == 2

    def test_correct_ratio(self):
        """Should compute gas/star ratio correctly."""
        result = ApertureFluxResult(
            star_flux=np.array([[10.0]]),
            gas_flux=np.array([[20.0]]),
            aperture_area_frac=np.ones((1, 1)),
            apertures=np.array([10.0]),
            peak_positions=np.array([[30, 30]]),
        )

        ratios, _ = compute_flux_ratios(result)

        assert np.isclose(ratios[0], 2.0)  # 20/10 = 2


class TestSelectNonOverlappingPeaks:
    """Tests for select_non_overlapping_peaks function."""

    def test_far_peaks_all_selected(self):
        """Well-separated peaks should all be selected."""
        peaks = np.array([
            [0, 0],
            [100, 0],
            [0, 100],
            [100, 100],
        ])

        selected = select_non_overlapping_peaks(peaks, min_separation=10.0)

        assert np.sum(selected) == 4

    def test_close_peaks_filtered(self):
        """Overlapping peaks should be filtered."""
        peaks = np.array([
            [0, 0],
            [5, 5],  # Close to first
            [100, 100],  # Far from both
        ])

        selected = select_non_overlapping_peaks(peaks, min_separation=20.0)

        # Either first or second should be selected, but not both
        assert np.sum(selected) <= 2

    def test_empty_input(self):
        """Empty input should return empty result."""
        peaks = np.array([]).reshape(0, 2)
        selected = select_non_overlapping_peaks(peaks, min_separation=10.0)
        assert len(selected) == 0

    def test_reproducibility_with_seed(self):
        """Same seed should give same result."""
        peaks = np.array([
            [0, 0],
            [5, 5],
            [10, 10],
        ])

        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(42)

        selected1 = select_non_overlapping_peaks(peaks, min_separation=8.0, rng=rng1)
        selected2 = select_non_overlapping_peaks(peaks, min_separation=8.0, rng=rng2)

        assert np.array_equal(selected1, selected2)


class TestMonteCarloFluxRatios:
    """Tests for monte_carlo_flux_ratios function."""

    def test_returns_correct_shapes(self):
        """Should return arrays of correct shape."""
        result = ApertureFluxResult(
            star_flux=np.array([[10.0, 20.0], [15.0, 25.0]]),
            gas_flux=np.array([[20.0, 40.0], [30.0, 50.0]]),
            aperture_area_frac=np.ones((2, 2)),
            apertures=np.array([10.0, 20.0]),
            peak_positions=np.array([[30, 30], [80, 80]]),  # Far apart
        )

        mean_ratios, std_ratios, n_peaks = monte_carlo_flux_ratios(
            result, n_mc=10, seed=42
        )

        assert len(mean_ratios) == 2
        assert len(std_ratios) == 2
        assert len(n_peaks) == 2

    def test_mean_is_reasonable(self):
        """Mean ratio should be close to simple ratio for well-separated peaks."""
        result = ApertureFluxResult(
            star_flux=np.array([[10.0, 20.0]]),
            gas_flux=np.array([[20.0, 40.0]]),
            aperture_area_frac=np.ones((1, 2)),
            apertures=np.array([5.0]),  # Small aperture
            peak_positions=np.array([[0, 0], [100, 100]]),  # Far apart
        )

        mean_ratios, _, _ = monte_carlo_flux_ratios(result, n_mc=50, seed=42)

        # Both peaks should always be selected, ratio should be (20+40)/(10+20) = 2
        assert np.isclose(mean_ratios[0], 2.0, rtol=0.1)

    def test_reproducibility(self):
        """Same seed should give same results."""
        result = ApertureFluxResult(
            star_flux=np.array([[10.0, 20.0, 30.0]]),
            gas_flux=np.array([[20.0, 40.0, 60.0]]),
            aperture_area_frac=np.ones((1, 3)),
            apertures=np.array([20.0]),
            peak_positions=np.array([[0, 0], [30, 30], [100, 100]]),
        )

        mean1, _, _ = monte_carlo_flux_ratios(result, n_mc=50, seed=42)
        mean2, _, _ = monte_carlo_flux_ratios(result, n_mc=50, seed=42)

        assert np.allclose(mean1, mean2)

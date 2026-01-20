"""
Unit tests for nearest-neighbour statistics module.

Tests the calculation of median relative nearest-neighbour distances
for peak analysis.
"""

import numpy as np
import pytest

from heisenberg.peaks.nearest_neighbour import (
    NearestNeighbourResult,
    med_peak_relative_nearest_neighbour_dist,
    nearest_neighbour_from_peaks,
)


class TestNearestNeighbourResult:
    """Tests for NearestNeighbourResult dataclass."""

    def test_create_result(self):
        """Should create result with all fields."""
        result = NearestNeighbourResult(
            dist_med=2.5,
            dist_med_sigma=0.3,
            dist_val=25.0,
            dist_sigma=3.0,
            fwhm_val=10.0,
            fwhm_sigma=1.0,
            n_peaks=5,
        )

        assert result.dist_med == 2.5
        assert result.dist_med_sigma == 0.3
        assert result.dist_val == 25.0
        assert result.dist_sigma == 3.0
        assert result.fwhm_val == 10.0
        assert result.fwhm_sigma == 1.0
        assert result.n_peaks == 5


class TestMedPeakRelativeNearestNeighbourDist:
    """Tests for med_peak_relative_nearest_neighbour_dist function."""

    def test_uniform_grid(self):
        """Uniform grid should have predictable NN distance."""
        # Create 3x3 grid of peaks at spacing 10
        x = np.array([0, 10, 20, 0, 10, 20, 0, 10, 20], dtype=float)
        y = np.array([0, 0, 0, 10, 10, 10, 20, 20, 20], dtype=float)
        fwhm_x = np.ones(9) * 5.0
        fwhm_y = np.ones(9) * 5.0

        result = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, seed=42
        )

        # NN distance should be 10 (grid spacing)
        assert result.dist_val == pytest.approx(10.0, rel=0.01)
        # Mean FWHM is 5.0
        assert result.fwhm_val == pytest.approx(5.0, rel=0.01)
        # Relative distance should be 10/5 = 2.0
        assert result.dist_med == pytest.approx(2.0, rel=0.01)
        assert result.n_peaks == 9

    def test_two_peaks(self):
        """Two peaks should have distance equal to their separation."""
        x = np.array([0.0, 30.0])
        y = np.array([0.0, 0.0])
        fwhm_x = np.array([5.0, 5.0])
        fwhm_y = np.array([5.0, 5.0])

        result = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, seed=42
        )

        assert result.dist_val == pytest.approx(30.0, rel=0.01)
        assert result.n_peaks == 2

    def test_single_peak_returns_nan(self):
        """Single peak should return NaN values."""
        x = np.array([10.0])
        y = np.array([10.0])
        fwhm_x = np.array([5.0])
        fwhm_y = np.array([5.0])

        result = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, seed=42
        )

        assert np.isnan(result.dist_med)
        assert np.isnan(result.dist_val)
        assert result.n_peaks == 1

    def test_no_peaks_returns_nan(self):
        """No peaks should return NaN values with n_peaks=0."""
        x = np.array([])
        y = np.array([])
        fwhm_x = np.array([])
        fwhm_y = np.array([])

        result = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, seed=42
        )

        assert np.isnan(result.dist_med)
        assert result.n_peaks == 0

    def test_reproducible_with_seed(self):
        """Same seed should produce same results."""
        x = np.array([0, 10, 25, 40, 50], dtype=float)
        y = np.array([0, 5, 15, 30, 45], dtype=float)
        fwhm_x = np.ones(5) * 5.0
        fwhm_y = np.ones(5) * 5.0

        result1 = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, seed=42
        )
        result2 = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, seed=42
        )

        assert result1.dist_med == result2.dist_med
        assert result1.dist_med_sigma == result2.dist_med_sigma

    def test_different_seeds_produce_similar_values(self):
        """Different seeds should produce similar median values."""
        x = np.array([0, 10, 25, 40, 50], dtype=float)
        y = np.array([0, 5, 15, 30, 45], dtype=float)
        fwhm_x = np.ones(5) * 5.0
        fwhm_y = np.ones(5) * 5.0

        result1 = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, seed=42
        )
        result2 = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, seed=123
        )

        # Median distance values should be identical (not affected by bootstrap)
        assert result1.dist_val == result2.dist_val

    def test_with_mask(self):
        """Mask should filter out peaks."""
        x = np.array([10, 50, 90], dtype=float)
        y = np.array([10, 50, 90], dtype=float)
        fwhm_x = np.ones(3) * 5.0
        fwhm_y = np.ones(3) * 5.0

        # Mask that excludes the middle peak
        mask = np.ones((100, 100))
        mask[45:55, 45:55] = 0

        result = med_peak_relative_nearest_neighbour_dist(
            x, y, fwhm_x, fwhm_y, mask=mask, seed=42
        )

        # Only 2 peaks should be used
        assert result.n_peaks == 2


class TestNearestNeighbourFromPeaks:
    """Tests for nearest_neighbour_from_peaks convenience function."""

    def test_with_mock_peaks(self):
        """Should work with DetectedPeak-like objects."""
        # Create mock peaks with stats
        class MockStats:
            def __init__(self, fwhm_x, fwhm_y):
                self.fwhm_x = fwhm_x
                self.fwhm_y = fwhm_y

        class MockPeak:
            def __init__(self, x, y, fwhm_x=5.0, fwhm_y=5.0):
                self.x = x
                self.y = y
                self.stats = MockStats(fwhm_x, fwhm_y)

        peaks = [
            MockPeak(0, 0),
            MockPeak(20, 0),
            MockPeak(40, 0),
        ]

        result = nearest_neighbour_from_peaks(peaks, seed=42)

        assert result.n_peaks == 3
        assert result.dist_val == pytest.approx(20.0, rel=0.01)

    def test_empty_peaks_list(self):
        """Empty peaks list should return NaN result."""
        result = nearest_neighbour_from_peaks([], seed=42)

        assert np.isnan(result.dist_med)
        assert result.n_peaks == 0

    def test_peaks_without_stats(self):
        """Peaks without stats should use default FWHM of 1.0."""
        class MockPeak:
            def __init__(self, x, y):
                self.x = x
                self.y = y
                self.stats = None

        peaks = [MockPeak(0, 0), MockPeak(10, 0), MockPeak(20, 0)]

        result = nearest_neighbour_from_peaks(peaks, seed=42)

        assert result.n_peaks == 3
        # FWHM should be 1.0 (default)
        assert result.fwhm_val == pytest.approx(1.0, rel=0.01)


class TestDistStatFromPeaks:
    """Tests for calculate_dist_stat_from_peaks in diffuse module."""

    def test_calculates_dist_stat(self):
        """Should calculate dist_stat normalized by cut_length."""
        from heisenberg.fourier.diffuse import calculate_dist_stat_from_peaks

        class MockStats:
            fwhm_x = 5.0
            fwhm_y = 5.0

        class MockPeak:
            def __init__(self, x, y):
                self.x = x
                self.y = y
                self.stats = MockStats()

        peaks = [MockPeak(0, 0), MockPeak(20, 0), MockPeak(40, 0)]
        cut_length = 10.0

        dist_stat, dist_stat_sigma = calculate_dist_stat_from_peaks(
            peaks, cut_length, seed=42
        )

        # NN distance is 20, cut_length is 10, so dist_stat = 2.0
        assert dist_stat == pytest.approx(2.0, rel=0.01)
        assert dist_stat_sigma >= 0

    def test_zero_cut_length_returns_nan(self):
        """Zero cut_length should return NaN."""
        from heisenberg.fourier.diffuse import calculate_dist_stat_from_peaks

        class MockPeak:
            def __init__(self, x, y):
                self.x = x
                self.y = y
                self.stats = None

        peaks = [MockPeak(0, 0), MockPeak(10, 0)]

        dist_stat, dist_stat_sigma = calculate_dist_stat_from_peaks(
            peaks, cut_length=0.0, seed=42
        )

        assert np.isnan(dist_stat)

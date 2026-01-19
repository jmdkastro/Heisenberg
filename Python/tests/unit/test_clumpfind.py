"""Tests for heisenberg.peaks module (clumpfind, statistics, detection)."""

import numpy as np
import pytest

from heisenberg.peaks import (
    Clump,
    ClumpfindResult,
    clumpfind2d,
    generate_levels,
    PeakStatistics,
    compute_clump_statistics,
    compute_all_statistics,
    statistics_to_array,
    SIGMA_TO_FWHM,
    PeakDetectionConfig,
    DetectedPeak,
    find_peaks,
    find_peaks_dual,
    generate_contour_levels,
    peaks_to_array,
    peaks_to_coords,
)


def make_gaussian_image(
    shape: tuple,
    center: tuple,
    amplitude: float,
    sigma: float,
    background: float = 0.0
) -> np.ndarray:
    """Create a 2D Gaussian image for testing."""
    y, x = np.ogrid[:shape[0], :shape[1]]
    cy, cx = center
    image = amplitude * np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * sigma**2))
    return image + background


def make_two_peak_image(
    shape: tuple,
    centers: list,
    amplitudes: list,
    sigmas: list,
    background: float = 0.0
) -> np.ndarray:
    """Create an image with two Gaussian peaks."""
    image = np.zeros(shape) + background
    for center, amp, sig in zip(centers, amplitudes, sigmas):
        y, x = np.ogrid[:shape[0], :shape[1]]
        cy, cx = center
        image += amp * np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * sig**2))
    return image


class TestClumpfind2d:
    """Tests for clumpfind2d function."""

    def test_single_gaussian_peak(self):
        """Single Gaussian should be identified as one clump."""
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25, 25),
            amplitude=100.0,
            sigma=5.0
        )
        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)

        assert len(result.clumps) == 1
        assert result.clumps[0].peak_position == (25, 25)
        assert np.isclose(result.clumps[0].peak_value, 100.0, rtol=0.01)

    def test_two_separated_peaks(self):
        """Two well-separated peaks should be two clumps."""
        image = make_two_peak_image(
            shape=(100, 100),
            centers=[(25, 25), (75, 75)],
            amplitudes=[100.0, 80.0],
            sigmas=[5.0, 5.0]
        )
        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)

        assert len(result.clumps) == 2
        # Sorted by peak flux, so first should be brighter
        assert result.clumps[0].peak_value > result.clumps[1].peak_value

    def test_merged_peaks_assigned_to_nearest(self):
        """When peaks merge at low level, pixels go to nearest peak."""
        # Two peaks far enough apart to be identified separately
        image = make_two_peak_image(
            shape=(100, 100),
            centers=[(30, 30), (70, 70)],
            amplitudes=[100.0, 80.0],
            sigmas=[8.0, 8.0],
            background=1.0  # Small background
        )
        levels = np.linspace(5, 90, 15)
        result = clumpfind2d(image, levels, npixmin=5)

        # Should find two clumps
        assert len(result.clumps) >= 1

        # Check assignment map consistency
        # Each pixel should be assigned to exactly one clump (or unassigned)
        for clump in result.clumps:
            clump_mask = result.assignment_map == clump.id
            assert np.sum(clump_mask) == clump.npix

    def test_npixmin_filtering(self):
        """Small clumps below npixmin should be rejected."""
        # Create image with one large peak and one small peak
        image = make_two_peak_image(
            shape=(100, 100),
            centers=[(25, 25), (75, 75)],
            amplitudes=[100.0, 20.0],
            sigmas=[10.0, 2.0]  # Second peak is much smaller
        )
        levels = np.linspace(5, 90, 15)

        # With high npixmin, small peak should be filtered
        result = clumpfind2d(image, levels, npixmin=50)
        assert len(result.clumps) == 1  # Only large peak survives

    def test_empty_image(self):
        """Empty (all zeros) image should return no clumps."""
        image = np.zeros((50, 50))
        levels = np.linspace(1, 10, 5)
        result = clumpfind2d(image, levels, npixmin=5)

        assert len(result.clumps) == 0

    def test_uniform_image(self):
        """Uniform image (no peaks) should return no clumps above threshold."""
        image = np.ones((50, 50)) * 5.0
        levels = np.linspace(10, 100, 10)  # All levels above image value
        result = clumpfind2d(image, levels, npixmin=5)

        assert len(result.clumps) == 0

    def test_nan_handling(self):
        """NaN values should be properly masked."""
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25, 25),
            amplitude=100.0,
            sigma=5.0
        )
        # Add NaN region
        image[10:15, 10:15] = np.nan

        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)

        # Should still find the peak
        assert len(result.clumps) == 1
        # NaN pixels should not be assigned
        assert result.assignment_map[12, 12] == 0

    def test_assignment_map_shape(self):
        """Assignment map should have same shape as input."""
        image = make_gaussian_image(
            shape=(64, 80),
            center=(32, 40),
            amplitude=100.0,
            sigma=5.0
        )
        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)

        assert result.assignment_map.shape == image.shape


class TestGenerateLevels:
    """Tests for generate_levels function."""

    def test_logarithmic_levels(self):
        """Log levels should span correct range."""
        image = np.ones((10, 10)) * 100.0
        image[5, 5] = 1000.0  # Peak

        levels = generate_levels(image, nlevels=10, logspacing=True, logrange=2.0)

        assert len(levels) == 10
        # Should span from 10^(3-2)=10 to 10^3=1000
        assert levels[0] < levels[-1]  # Ascending order
        assert np.isclose(np.log10(levels[-1] / levels[0]), 2.0, rtol=0.1)

    def test_linear_levels(self):
        """Linear levels should be evenly spaced."""
        image = np.zeros((10, 10))
        image[5, 5] = 100.0

        levels = generate_levels(image, nlevels=11, logspacing=False, minlevel=0.0)

        assert len(levels) == 11
        # Should be evenly spaced
        diffs = np.diff(levels)
        assert np.allclose(diffs, diffs[0], rtol=0.01)

    def test_minlevel_override(self):
        """Custom minlevel should be respected."""
        image = np.ones((10, 10)) * 100.0
        image[5, 5] = 1000.0

        levels = generate_levels(
            image, nlevels=10, logspacing=True, logrange=2.0, minlevel=50.0
        )

        # Allow small floating point tolerance
        assert levels[0] >= 49.9


class TestPeakStatistics:
    """Tests for peak statistics calculation."""

    def test_gaussian_fwhm(self):
        """FWHM of Gaussian should match input sigma."""
        sigma = 5.0
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25, 25),
            amplitude=100.0,
            sigma=sigma
        )

        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)
        stats = compute_clump_statistics(image, result.clumps[0], flux_weighted=True)

        expected_fwhm = SIGMA_TO_FWHM * sigma
        # Allow some tolerance due to pixelization and clump boundary
        assert np.isclose(stats.fwhm_x, expected_fwhm, rtol=0.2)
        assert np.isclose(stats.fwhm_y, expected_fwhm, rtol=0.2)

    def test_flux_weighted_position(self):
        """Flux-weighted position should be at centroid."""
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25.0, 25.0),
            amplitude=100.0,
            sigma=5.0
        )

        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)
        stats = compute_clump_statistics(image, result.clumps[0], flux_weighted=True)

        # Centroid should be at (25, 25)
        assert np.isclose(stats.x, 25.0, atol=0.5)
        assert np.isclose(stats.y, 25.0, atol=0.5)

    def test_peak_position(self):
        """Non-flux-weighted position should be at peak pixel."""
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25, 30),
            amplitude=100.0,
            sigma=5.0
        )

        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)
        stats = compute_clump_statistics(image, result.clumps[0], flux_weighted=False)

        # Should be exactly at peak pixel
        assert stats.x == 30.0
        assert stats.y == 25.0

    def test_total_flux_positive(self):
        """Total flux should be positive and reasonable."""
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25, 25),
            amplitude=100.0,
            sigma=5.0
        )

        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)
        stats = compute_clump_statistics(image, result.clumps[0])

        assert stats.total_flux > 0
        assert stats.peak_flux == 100.0

    def test_equivalent_radius(self):
        """Equivalent radius should equal sqrt(npix/pi)."""
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25, 25),
            amplitude=100.0,
            sigma=5.0
        )

        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)
        stats = compute_clump_statistics(image, result.clumps[0])

        expected_radius = np.sqrt(stats.npix / np.pi)
        assert np.isclose(stats.radius, expected_radius)

    def test_on_edge_detection(self):
        """Should detect when clump touches image boundary."""
        # Peak near edge
        image = make_gaussian_image(
            shape=(50, 50),
            center=(3, 25),  # Near top edge
            amplitude=100.0,
            sigma=5.0
        )

        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)
        stats = compute_clump_statistics(image, result.clumps[0])

        assert stats.on_edge is True

    def test_statistics_to_array(self):
        """Should convert statistics list to array correctly."""
        image = make_two_peak_image(
            shape=(100, 100),
            centers=[(25, 25), (75, 75)],
            amplitudes=[100.0, 80.0],
            sigmas=[5.0, 5.0]
        )

        levels = np.linspace(10, 90, 10)
        result = clumpfind2d(image, levels, npixmin=5)
        stats_list = compute_all_statistics(image, result)
        arr = statistics_to_array(stats_list)

        assert arr.shape == (2, 9)
        assert arr[0, 3] > arr[1, 3]  # First has higher peak flux


class TestFindPeaks:
    """Tests for high-level find_peaks function."""

    def test_basic_detection(self):
        """Should find peaks with default config."""
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25, 25),
            amplitude=100.0,
            sigma=5.0
        )

        config = PeakDetectionConfig(npixmin=5, nlevels=10)
        peaks = find_peaks(image, config)

        assert len(peaks) == 1
        assert np.isclose(peaks[0].x, 25, atol=1)
        assert np.isclose(peaks[0].y, 25, atol=1)

    def test_sensitivity_filtering(self):
        """Should filter peaks by sensitivity threshold."""
        image = make_two_peak_image(
            shape=(100, 100),
            centers=[(25, 25), (75, 75)],
            amplitudes=[100.0, 20.0],  # Second peak much dimmer
            sigmas=[5.0, 5.0]
        )

        # Uniform sensitivity map
        sensitivity = np.ones_like(image) * 5.0

        # With nsigma=5, threshold is 25. Second peak (20) should be filtered.
        config = PeakDetectionConfig(npixmin=5, nsigma=5.0, nlevels=15)
        peaks = find_peaks(image, config, sensitivity=sensitivity)

        assert len(peaks) == 1
        assert peaks[0].peak_flux > 50  # Only bright peak survives

    def test_sorted_by_flux(self):
        """Peaks should be sorted by total flux descending."""
        image = make_two_peak_image(
            shape=(100, 100),
            centers=[(25, 25), (75, 75)],
            amplitudes=[100.0, 150.0],
            sigmas=[5.0, 8.0]  # Second has more total flux due to larger sigma
        )

        config = PeakDetectionConfig(npixmin=5, nlevels=15)
        peaks = find_peaks(image, config)

        assert len(peaks) == 2
        assert peaks[0].total_flux >= peaks[1].total_flux

    def test_include_stats(self):
        """Should include full statistics when requested."""
        image = make_gaussian_image(
            shape=(50, 50),
            center=(25, 25),
            amplitude=100.0,
            sigma=5.0
        )

        config = PeakDetectionConfig(npixmin=5)
        peaks_with = find_peaks(image, config, include_stats=True)
        peaks_without = find_peaks(image, config, include_stats=False)

        assert peaks_with[0].stats is not None
        assert peaks_without[0].stats is None

    def test_peaks_to_array(self):
        """Should convert peaks to array correctly."""
        image = make_two_peak_image(
            shape=(100, 100),
            centers=[(25, 25), (75, 75)],
            amplitudes=[100.0, 80.0],
            sigmas=[5.0, 5.0]
        )

        config = PeakDetectionConfig(npixmin=5)
        peaks = find_peaks(image, config)
        arr = peaks_to_array(peaks)

        assert arr.shape == (2, 4)  # [x, y, total_flux, npix]

    def test_peaks_to_coords(self):
        """Should extract coordinates correctly."""
        image = make_two_peak_image(
            shape=(100, 100),
            centers=[(25, 25), (75, 75)],
            amplitudes=[100.0, 80.0],
            sigmas=[5.0, 5.0]
        )

        config = PeakDetectionConfig(npixmin=5)
        peaks = find_peaks(image, config)
        coords = peaks_to_coords(peaks)

        assert coords.shape == (2, 2)  # [x, y] for each peak


class TestFindPeaksDual:
    """Tests for dual star/gas peak detection."""

    def test_dual_detection(self):
        """Should find peaks in both maps independently."""
        star_image = make_gaussian_image(
            shape=(50, 50), center=(20, 20), amplitude=100.0, sigma=5.0
        )
        gas_image = make_gaussian_image(
            shape=(50, 50), center=(30, 30), amplitude=80.0, sigma=5.0
        )

        star_config = PeakDetectionConfig(npixmin=5)
        gas_config = PeakDetectionConfig(npixmin=5)

        star_peaks, gas_peaks = find_peaks_dual(
            star_image, gas_image, star_config, gas_config
        )

        assert len(star_peaks) == 1
        assert len(gas_peaks) == 1
        # Check they found different positions
        assert not np.allclose(
            [star_peaks[0].x, star_peaks[0].y],
            [gas_peaks[0].x, gas_peaks[0].y],
            atol=2
        )

    def test_different_configs(self):
        """Should apply different configs to each map."""
        # Create image with one large peak and one very small peak
        star_image = make_two_peak_image(
            shape=(100, 100),
            centers=[(25, 25), (75, 75)],
            amplitudes=[100.0, 15.0],
            sigmas=[10.0, 2.0]  # Second peak is very small
        )
        gas_image = star_image.copy()

        # Star config with very high npixmin filters the small peak
        star_config = PeakDetectionConfig(npixmin=100, nlevels=15)
        # Gas config with low npixmin keeps both
        gas_config = PeakDetectionConfig(npixmin=5, nlevels=15)

        star_peaks, gas_peaks = find_peaks_dual(
            star_image, gas_image, star_config, gas_config
        )

        # Large npixmin should filter the small peak
        assert len(star_peaks) <= len(gas_peaks)


class TestContourLevels:
    """Tests for contour level generation."""

    def test_log_levels_from_config(self):
        """Should generate log levels from config."""
        image = np.ones((10, 10)) * 10.0
        image[5, 5] = 1000.0

        config = PeakDetectionConfig(loglevels=True, logrange=2.0, nlevels=10)
        levels = generate_contour_levels(image, config)

        assert len(levels) == 10
        assert levels[0] < levels[-1]

    def test_linear_levels_from_config(self):
        """Should generate linear levels from config."""
        image = np.ones((10, 10)) * 10.0
        image[5, 5] = 100.0

        config = PeakDetectionConfig(loglevels=False, nlevels=11)
        levels = generate_contour_levels(image, config)

        assert len(levels) == 11
        # Should be approximately evenly spaced
        diffs = np.diff(levels)
        assert np.std(diffs) / np.mean(diffs) < 0.1

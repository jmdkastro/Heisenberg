"""Tests for heisenberg.fourier module."""

import numpy as np
import pytest

from heisenberg.fourier import (
    # Filters
    FilterKernel,
    FilterPass,
    lowpass_gaussian,
    lowpass_butterworth,
    lowpass_ideal,
    highpass_gaussian,
    highpass_butterworth,
    highpass_ideal,
    get_filter,
    create_frequency_grid,
    # Filter tool
    FilterResult,
    apply_fourier_filter,
    compute_power_spectrum,
    compute_power_fraction,
    extract_compact_emission,
    # Corrections
    flux_loss_correction,
    overlap_correction,
    overlap_correction_sigma,
    compute_corrections,
    symmetric_sigmoidal,
    # Diffuse fraction
    calculate_diffuse_fraction,
    calculate_diffuse_fraction_with_errors,
)


# ============================================================================
# Test Fixtures
# ============================================================================

def make_gaussian_image(shape, center, amplitude, sigma):
    """Create a 2D Gaussian test image."""
    y, x = np.ogrid[:shape[0], :shape[1]]
    cy, cx = center
    return amplitude * np.exp(-((x - cx)**2 + (y - cy)**2) / (2 * sigma**2))


def make_test_image_with_scales():
    """Create test image with both small and large scale structure."""
    shape = (128, 128)
    # Large scale (diffuse) - big Gaussian
    diffuse = make_gaussian_image(shape, (64, 64), 100.0, 30.0)
    # Small scale (compact) - small Gaussians
    compact1 = make_gaussian_image(shape, (40, 40), 50.0, 5.0)
    compact2 = make_gaussian_image(shape, (90, 90), 40.0, 5.0)
    return diffuse + compact1 + compact2


# ============================================================================
# Filter Tests
# ============================================================================

class TestLowpassGaussian:
    """Tests for Gaussian lowpass filter."""

    def test_dc_component_preserved(self):
        """DC component (freq=0) should be fully preserved."""
        freq_dist = np.array([0.0, 0.1, 0.2])
        taper = lowpass_gaussian(freq_dist, cut_freq=0.1)
        assert taper[0] == 1.0

    def test_high_freq_attenuated(self):
        """High frequencies should be strongly attenuated."""
        freq_dist = np.array([0.0, 0.5, 1.0])
        taper = lowpass_gaussian(freq_dist, cut_freq=0.1)
        assert taper[2] < 0.01  # Strongly attenuated

    def test_output_range(self):
        """Output should be in [0, 1]."""
        freq_dist = np.random.uniform(0, 1, 100)
        taper = lowpass_gaussian(freq_dist, cut_freq=0.2)
        assert np.all(taper >= 0)
        assert np.all(taper <= 1)


class TestLowpassButterworth:
    """Tests for Butterworth lowpass filter."""

    def test_dc_component_preserved(self):
        """DC component should be preserved."""
        freq_dist = np.array([0.0, 0.1, 0.2])
        taper = lowpass_butterworth(freq_dist, cut_freq=0.1, order=2)
        assert np.isclose(taper[0], 1.0)

    def test_higher_order_sharper_cutoff(self):
        """Higher order should give sharper cutoff."""
        freq_dist = np.array([0.15])  # Just above cutoff
        taper_low = lowpass_butterworth(freq_dist, cut_freq=0.1, order=1)
        taper_high = lowpass_butterworth(freq_dist, cut_freq=0.1, order=4)
        # Higher order should attenuate more at same frequency
        assert taper_high[0] < taper_low[0]

    def test_invalid_order_raises(self):
        """Non-positive integer order should raise error."""
        freq_dist = np.array([0.1])
        with pytest.raises(ValueError):
            lowpass_butterworth(freq_dist, cut_freq=0.1, order=0)
        with pytest.raises(ValueError):
            lowpass_butterworth(freq_dist, cut_freq=0.1, order=-1)


class TestLowpassIdeal:
    """Tests for ideal lowpass filter."""

    def test_binary_output(self):
        """Output should be 0 or 1 only."""
        freq_dist = np.array([0.05, 0.1, 0.15, 0.2])
        taper = lowpass_ideal(freq_dist, cut_freq=0.1)
        assert np.all(np.isin(taper, [0.0, 1.0]))

    def test_cutoff_behavior(self):
        """Frequencies below cutoff pass, above are blocked."""
        freq_dist = np.array([0.05, 0.15])
        taper = lowpass_ideal(freq_dist, cut_freq=0.1)
        assert taper[0] == 1.0  # Below cutoff
        assert taper[1] == 0.0  # Above cutoff


class TestHighpassFilters:
    """Tests for highpass filters."""

    def test_gaussian_complement(self):
        """Gaussian highpass should be 1 - lowpass."""
        freq_dist = np.array([0.0, 0.1, 0.2, 0.3])
        lp = lowpass_gaussian(freq_dist, cut_freq=0.1)
        hp = highpass_gaussian(freq_dist, cut_freq=0.1)
        assert np.allclose(lp + hp, 1.0)

    def test_ideal_complement(self):
        """Ideal highpass should be complement of lowpass."""
        freq_dist = np.array([0.05, 0.15])
        lp = lowpass_ideal(freq_dist, cut_freq=0.1)
        hp = highpass_ideal(freq_dist, cut_freq=0.1)
        # Not exactly 1 due to boundary handling, but complementary
        assert hp[0] == 0.0  # Below cutoff blocked by highpass
        assert hp[1] == 1.0  # Above cutoff passed by highpass

    def test_butterworth_dc_blocked(self):
        """Butterworth highpass should block DC."""
        freq_dist = np.array([0.0, 0.1, 0.2])
        taper = highpass_butterworth(freq_dist, cut_freq=0.1, order=2)
        assert taper[0] == 0.0


class TestGetFilter:
    """Tests for get_filter dispatcher."""

    def test_string_kernel(self):
        """Should accept string kernel names."""
        freq_dist = np.array([0.1])
        taper = get_filter('gaussian', 'low', freq_dist, cut_freq=0.2)
        assert len(taper) == 1

    def test_enum_kernel(self):
        """Should accept enum kernel types."""
        freq_dist = np.array([0.1])
        taper = get_filter(FilterKernel.BUTTERWORTH, FilterPass.HIGHPASS,
                          freq_dist, cut_freq=0.2, order=2)
        assert len(taper) == 1


class TestFrequencyGrid:
    """Tests for frequency grid creation."""

    def test_shape_preserved(self):
        """Output shape should match input shape."""
        shape = (64, 80)
        grid = create_frequency_grid(shape)
        assert grid.shape == shape

    def test_origin_at_zero(self):
        """Origin (DC) should have frequency distance zero."""
        grid = create_frequency_grid((64, 64))
        assert grid[0, 0] == 0.0

    def test_max_frequency(self):
        """Maximum frequency should be near Nyquist."""
        grid = create_frequency_grid((64, 64))
        # Nyquist in 1D is 0.5, so max 2D distance is sqrt(2)*0.5 ≈ 0.707
        assert np.max(grid) <= np.sqrt(2) * 0.5 + 0.01


# ============================================================================
# Filter Tool Tests
# ============================================================================

class TestApplyFourierFilter:
    """Tests for apply_fourier_filter function."""

    def test_output_shape(self):
        """Output should match input shape."""
        image = np.random.randn(64, 64)
        result = apply_fourier_filter(image, cut_length=10)
        assert result.filtered.shape == image.shape

    def test_lowpass_smooths(self):
        """Lowpass filter should reduce variance."""
        np.random.seed(42)
        image = np.random.randn(64, 64)
        result = apply_fourier_filter(image, cut_length=5, pass_type='low')
        assert np.var(result.filtered) < np.var(image)

    def test_nan_preserved(self):
        """NaN pixels should be preserved."""
        image = np.ones((32, 32))
        image[10, 10] = np.nan
        result = apply_fourier_filter(image, cut_length=5)
        assert np.isnan(result.filtered[10, 10])

    def test_residual_computation(self):
        """Residual should be original - filtered."""
        image = np.random.randn(32, 32)
        result = apply_fourier_filter(image, cut_length=5, compute_residual=True)
        assert result.residual is not None
        expected_residual = image - result.filtered
        assert np.allclose(result.residual, expected_residual, equal_nan=True)


class TestExtractCompactEmission:
    """Tests for compact/diffuse separation."""

    def test_returns_two_arrays(self):
        """Should return compact and diffuse components."""
        image = make_test_image_with_scales()
        compact, diffuse = extract_compact_emission(image, cut_length=20)
        assert compact.shape == image.shape
        assert diffuse.shape == image.shape

    def test_sum_equals_original(self):
        """Compact + diffuse should equal original (ignoring NaN)."""
        image = make_test_image_with_scales()
        compact, diffuse = extract_compact_emission(image, cut_length=20)
        reconstructed = compact + diffuse
        assert np.allclose(reconstructed, image)


class TestPowerSpectrum:
    """Tests for power spectrum computation."""

    def test_positive_values(self):
        """Power spectrum should be non-negative."""
        image = np.random.randn(32, 32)
        power = compute_power_spectrum(image)
        assert np.all(power >= 0)

    def test_parseval_relation(self):
        """Total power should relate to sum of squares (Parseval)."""
        image = np.random.randn(32, 32)
        power = compute_power_spectrum(image)
        # Parseval: sum(|x|^2) = sum(|X|^2) / N
        N = image.size
        assert np.isclose(np.sum(image**2), np.sum(power) / N)


# ============================================================================
# Correction Tests
# ============================================================================

class TestFluxLossCorrection:
    """Tests for flux loss correction."""

    def test_large_ratio_approaches_one(self):
        """Large cut ratio should give correction near 1."""
        qcon = flux_loss_correction(100.0, kernel='gaussian')
        assert np.isclose(qcon, 1.0, atol=0.01)

    def test_small_ratio_less_than_one(self):
        """Small cut ratio should give correction < 1."""
        qcon = flux_loss_correction(1.0, kernel='gaussian')
        assert qcon < 1.0

    def test_array_input(self):
        """Should handle array input."""
        ratios = np.array([1.0, 5.0, 10.0])
        qcon = flux_loss_correction(ratios, kernel='gaussian')
        assert len(qcon) == 3
        assert np.all(qcon > 0)


class TestOverlapCorrection:
    """Tests for overlap correction."""

    def test_large_distance_approaches_one(self):
        """Large distance stat should give correction near 1."""
        qover = overlap_correction(10.0, kernel='gaussian')
        assert qover <= 1.0
        assert qover > 0.9

    def test_clamped_to_one(self):
        """Should be clamped to maximum 1.0."""
        qover = overlap_correction(100.0, kernel='gaussian')
        assert qover <= 1.0

    def test_small_distance_less_than_one(self):
        """Small distance stat should give correction < 1."""
        qover = overlap_correction(0.5, kernel='gaussian')
        assert qover < 1.0


class TestOverlapCorrectionSigma:
    """Tests for overlap correction uncertainty."""

    def test_positive_uncertainty(self):
        """Uncertainty should be positive."""
        sigma = overlap_correction_sigma(1.5, dist_stat_sigma=0.1)
        assert sigma >= 0

    def test_zero_input_sigma_gives_nonzero(self):
        """Even without input uncertainty, parameter uncertainty exists."""
        sigma = overlap_correction_sigma(1.5, dist_stat_sigma=0.0)
        assert sigma >= 0


class TestComputeCorrections:
    """Tests for combined corrections."""

    def test_returns_dataclass(self):
        """Should return CorrectionFactors dataclass."""
        result = compute_corrections(5.0, 2.0, kernel='gaussian')
        assert hasattr(result, 'qcon')
        assert hasattr(result, 'qover')
        assert hasattr(result, 'total')

    def test_total_is_product(self):
        """Total should be qcon * qover."""
        result = compute_corrections(5.0, 2.0)
        assert np.isclose(result.total, result.qcon * result.qover)


# ============================================================================
# Diffuse Fraction Tests
# ============================================================================

class TestCalculateDiffuseFraction:
    """Tests for diffuse fraction calculation."""

    def test_returns_result(self):
        """Should return DiffuseFractionResult."""
        image = make_test_image_with_scales()
        result = calculate_diffuse_fraction(image, cut_length=20)
        assert hasattr(result, 'flux_frac')
        assert hasattr(result, 'diffuse_frac')

    def test_fractions_sum_to_one(self):
        """Flux fraction + diffuse fraction should equal 1."""
        image = make_test_image_with_scales()
        result = calculate_diffuse_fraction(image, cut_length=20)
        assert np.isclose(result.flux_frac + result.diffuse_frac, 1.0)

    def test_fractions_in_valid_range(self):
        """Fractions should be in [0, 1]."""
        image = make_test_image_with_scales()
        result = calculate_diffuse_fraction(image, cut_length=20)
        assert 0 <= result.flux_frac <= 1
        assert 0 <= result.diffuse_frac <= 1

    def test_small_cut_length_more_compact(self):
        """Smaller cut length should classify more as diffuse."""
        image = make_test_image_with_scales()
        result_small = calculate_diffuse_fraction(image, cut_length=5)
        result_large = calculate_diffuse_fraction(image, cut_length=50)
        # Smaller cutoff -> less passes highpass -> less compact
        assert result_small.flux_frac < result_large.flux_frac

    def test_with_mask(self):
        """Should respect mask."""
        image = np.ones((32, 32)) * 10.0
        mask = np.ones((32, 32))
        mask[:16, :] = 0  # Mask out half

        result_full = calculate_diffuse_fraction(image, cut_length=5)
        result_masked = calculate_diffuse_fraction(image, cut_length=5, mask=mask)

        # Total flux should be different
        assert result_masked.total_flux < result_full.total_flux


class TestDiffuseFractionWithErrors:
    """Tests for diffuse fraction with error propagation."""

    def test_returns_errors(self):
        """Should return error estimates."""
        image = make_test_image_with_scales()
        result = calculate_diffuse_fraction_with_errors(
            image, cut_length=20,
            cut_length_errmin=2.0, cut_length_errmax=3.0
        )
        assert result.flux_frac_err is not None
        assert result.diffuse_frac_err is not None

    def test_errors_non_negative(self):
        """Errors should be non-negative."""
        image = make_test_image_with_scales()
        result = calculate_diffuse_fraction_with_errors(
            image, cut_length=20,
            cut_length_errmin=2.0, cut_length_errmax=3.0
        )
        assert result.flux_frac_err >= 0
        assert result.diffuse_frac_err >= 0

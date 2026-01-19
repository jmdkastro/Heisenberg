"""Tests for heisenberg.imaging module."""

import numpy as np
import pytest
from pathlib import Path
import tempfile

from heisenberg.imaging import (
    psf_gaussian,
    psf_tophat,
    convolve_image,
    Mask,
    create_circular_mask,
    create_radial_mask,
    mask_from_nan,
    synchronize_masks,
)


class TestPsfGaussian:
    """Tests for psf_gaussian function."""

    def test_normalized(self):
        """Gaussian kernel should sum to 1."""
        kernel = psf_gaussian(fwhm=5.0)
        assert np.isclose(np.sum(kernel), 1.0, rtol=1e-5)

    def test_peak_at_center(self):
        """Maximum should be at the center."""
        kernel = psf_gaussian(fwhm=5.0, size=11)
        center = kernel.shape[0] // 2
        assert kernel[center, center] == np.max(kernel)

    def test_symmetric(self):
        """Kernel should be symmetric."""
        kernel = psf_gaussian(fwhm=5.0, size=11)
        assert np.allclose(kernel, kernel.T)
        assert np.allclose(kernel, np.flip(kernel, axis=0))

    def test_size_odd(self):
        """Kernel size should be odd."""
        kernel = psf_gaussian(fwhm=5.0)
        assert kernel.shape[0] % 2 == 1
        assert kernel.shape[1] % 2 == 1


class TestPsfTophat:
    """Tests for psf_tophat function."""

    def test_normalized(self):
        """Tophat kernel should sum to 1."""
        kernel = psf_tophat(width=10.0)
        assert np.isclose(np.sum(kernel), 1.0, rtol=1e-3)

    def test_size_odd(self):
        """Kernel size should be odd."""
        kernel = psf_tophat(width=10.0)
        assert kernel.shape[0] % 2 == 1
        assert kernel.shape[1] % 2 == 1

    def test_center_nonzero(self):
        """Center pixel should be nonzero."""
        kernel = psf_tophat(width=5.0)
        center = kernel.shape[0] // 2
        assert kernel[center, center] > 0

    def test_corners_zero_for_large_width(self):
        """Corner pixels should be zero for circular kernel."""
        kernel = psf_tophat(width=10.0)
        # Corners should be outside the circle
        assert kernel[0, 0] == 0


class TestConvolveImage:
    """Tests for convolve_image function."""

    def test_preserves_shape(self):
        """Output should have same shape as input."""
        image = np.random.randn(50, 50)
        kernel = psf_gaussian(fwhm=3.0, size=7)
        result = convolve_image(image, kernel)
        assert result.shape == image.shape

    def test_preserves_nan(self):
        """NaN pixels should be preserved."""
        image = np.ones((20, 20))
        image[10, 10] = np.nan
        kernel = psf_gaussian(fwhm=2.0, size=5)
        result = convolve_image(image, kernel, preserve_nan=True)
        assert np.isnan(result[10, 10])

    def test_smoothing_reduces_variance(self):
        """Smoothing should reduce image variance."""
        np.random.seed(42)
        image = np.random.randn(50, 50)
        kernel = psf_gaussian(fwhm=5.0, size=15)
        result = convolve_image(image, kernel)
        assert np.var(result) < np.var(image)


class TestMask:
    """Tests for Mask class."""

    def test_apply_mask(self):
        """Applying mask should set masked pixels to NaN."""
        image = np.ones((10, 10))
        mask = Mask(data=np.array([[1, 0], [0, 1]]).repeat(5, axis=0).repeat(5, axis=1).astype(float))
        # Resize mask to match
        mask_data = np.ones((10, 10))
        mask_data[0:5, 5:10] = 0
        mask_data[5:10, 0:5] = 0
        mask = Mask(data=mask_data)
        result = mask.apply(image)
        assert np.isnan(result[2, 7])  # Masked region
        assert result[2, 2] == 1.0  # Unmasked region

    def test_combine_masks(self):
        """Combining masks should use logical AND."""
        mask1 = Mask(data=np.array([[1, 1], [0, 0]], dtype=float))
        mask2 = Mask(data=np.array([[1, 0], [1, 0]], dtype=float))
        combined = mask1.combine(mask2)
        expected = np.array([[1, 0], [0, 0]], dtype=float)
        assert np.array_equal(combined.data, expected)

    def test_n_unmasked(self):
        """n_unmasked should count unmasked pixels."""
        mask = Mask(data=np.array([[1, 1], [0, 1]], dtype=float))
        assert mask.n_unmasked == 3
        assert mask.n_masked == 1


class TestCreateCircularMask:
    """Tests for create_circular_mask function."""

    def test_mask_inside(self):
        """Mask inside circle should mask center."""
        mask = create_circular_mask(
            shape=(21, 21),
            center=(10, 10),
            radius=5,
            inside=True
        )
        # Center should be masked (0)
        assert mask.data[10, 10] == 0
        # Far corner should be unmasked (1)
        assert mask.data[0, 0] == 1

    def test_mask_outside(self):
        """Mask outside circle should mask corners."""
        mask = create_circular_mask(
            shape=(21, 21),
            center=(10, 10),
            radius=5,
            inside=False
        )
        # Center should be unmasked (1)
        assert mask.data[10, 10] == 1
        # Far corner should be masked (0)
        assert mask.data[0, 0] == 0


class TestCreateRadialMask:
    """Tests for create_radial_mask function."""

    def test_annulus_mask(self):
        """Radial mask should create an annulus."""
        mask = create_radial_mask(
            shape=(101, 101),
            center=(50, 50),
            inner_radius=10,
            outer_radius=30,
            pixel_to_pc=1.0
        )
        # Center should be masked (inside inner radius)
        assert mask.data[50, 50] == 0
        # Within annulus should be unmasked
        assert mask.data[50, 70] == 1  # 20 pixels from center
        # Outside outer radius should be masked
        assert mask.data[50, 90] == 0  # 40 pixels from center


class TestMaskFromNan:
    """Tests for mask_from_nan function."""

    def test_masks_nan_pixels(self):
        """NaN pixels should be masked."""
        image = np.ones((5, 5))
        image[2, 2] = np.nan
        image[0, 0] = np.nan
        mask = mask_from_nan(image)
        assert mask.data[2, 2] == 0
        assert mask.data[0, 0] == 0
        assert mask.data[1, 1] == 1


class TestSynchronizeMasks:
    """Tests for synchronize_masks function."""

    def test_combines_nan_masks(self):
        """Should mask pixels that are NaN in any image."""
        img1 = np.ones((5, 5))
        img1[0, 0] = np.nan
        img2 = np.ones((5, 5))
        img2[1, 1] = np.nan
        mask = synchronize_masks([img1, img2])
        assert mask.data[0, 0] == 0  # NaN in img1
        assert mask.data[1, 1] == 0  # NaN in img2
        assert mask.data[2, 2] == 1  # Valid in both

    def test_requires_same_shape(self):
        """Should raise error for different shapes."""
        img1 = np.ones((5, 5))
        img2 = np.ones((6, 6))
        with pytest.raises(ValueError):
            synchronize_masks([img1, img2])

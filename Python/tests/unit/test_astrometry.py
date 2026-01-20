"""
Unit tests for astrometry module.

Tests WCS comparison and astrometric utilities.
"""

import numpy as np
import pytest

from heisenberg.imaging.astrometry import (
    get_equinox,
    get_platescale,
    get_rotation,
    astrometry_equal,
    images_aligned,
    compute_pixel_offset,
)


class TestGetEquinox:
    """Tests for get_equinox function."""

    def test_equinox_keyword(self):
        """Should extract equinox from EQUINOX keyword."""
        header = {'EQUINOX': 2000.0}
        assert get_equinox(header) == 2000.0

    def test_epoch_keyword(self):
        """Should extract equinox from EPOCH keyword."""
        header = {'EPOCH': 1950.0}
        assert get_equinox(header) == 1950.0

    def test_radesys_fk5(self):
        """Should infer J2000 from FK5 RADESYS."""
        header = {'RADESYS': 'FK5'}
        assert get_equinox(header) == 2000.0

    def test_radesys_icrs(self):
        """Should infer J2000 from ICRS RADESYS."""
        header = {'RADESYS': 'ICRS'}
        assert get_equinox(header) == 2000.0

    def test_radesys_fk4(self):
        """Should infer B1950 from FK4 RADESYS."""
        header = {'RADESYS': 'FK4'}
        assert get_equinox(header) == 1950.0

    def test_no_equinox_info(self):
        """Should return None when no equinox info available."""
        header = {'NAXIS': 2}
        assert get_equinox(header) is None


class TestGetPlatescale:
    """Tests for get_platescale function."""

    def test_with_cdelt(self):
        """Should extract plate scale from CDELT."""
        from astropy.wcs import WCS

        # Create WCS with CDELT
        header = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,  # degrees per pixel
            'CDELT2': 0.001,
        }
        wcs = WCS(header)
        cdelt_x, cdelt_y = get_platescale(wcs)

        assert abs(cdelt_x) == pytest.approx(0.001, rel=0.01)
        assert abs(cdelt_y) == pytest.approx(0.001, rel=0.01)

    def test_raises_without_wcs(self):
        """Should raise when WCS is None."""
        with pytest.raises(ValueError, match="No WCS provided"):
            get_platescale(None)


class TestGetRotation:
    """Tests for get_rotation function."""

    def test_no_rotation(self):
        """Should return 0 for unrotated WCS."""
        from astropy.wcs import WCS

        header = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
        }
        wcs = WCS(header)
        rotation = get_rotation(wcs)

        # Should be close to 0 or 180 for standard orientation
        assert abs(rotation) < 1.0 or abs(abs(rotation) - 180) < 1.0

    def test_raises_without_wcs(self):
        """Should raise when WCS is None."""
        with pytest.raises(ValueError, match="No WCS provided"):
            get_rotation(None)


class TestAstrometryEqual:
    """Tests for astrometry_equal function."""

    def test_identical_images(self):
        """Identical images should have equal astrometry."""
        from astropy.wcs import WCS

        header = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
            'EQUINOX': 2000.0,
        }

        image = np.random.rand(100, 100)

        assert astrometry_equal(image, header, image, header) is True

    def test_different_shapes(self):
        """Different shaped images should not be equal."""
        header = {
            'NAXIS': 2,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
        }

        image_a = np.random.rand(100, 100)
        image_b = np.random.rand(50, 50)

        assert astrometry_equal(image_a, header, image_b, header) is False

    def test_different_crval(self):
        """Images with different CRVAL should not be equal."""
        from astropy.wcs import WCS

        header_a = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
        }

        header_b = header_a.copy()
        header_b['CRVAL1'] = 181.0  # Different RA

        image = np.random.rand(100, 100)

        assert astrometry_equal(image, header_a, image, header_b) is False

    def test_different_plate_scale(self):
        """Images with different plate scales should not be equal."""
        header_a = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
        }

        header_b = header_a.copy()
        header_b['CDELT1'] = -0.002  # Different plate scale

        image = np.random.rand(100, 100)

        assert astrometry_equal(image, header_a, image, header_b) is False


class TestImagesAligned:
    """Tests for images_aligned function."""

    def test_identical_wcs(self):
        """Identical WCS should be aligned."""
        from astropy.wcs import WCS

        header = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
        }

        wcs = WCS(header)
        image = np.random.rand(100, 100)

        assert images_aligned(image, wcs, image, wcs) is True

    def test_different_shapes_not_aligned(self):
        """Different shaped images should not be aligned."""
        from astropy.wcs import WCS

        header = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
        }

        wcs = WCS(header)
        image_a = np.random.rand(100, 100)
        image_b = np.random.rand(50, 50)

        assert images_aligned(image_a, wcs, image_b, wcs) is False


class TestComputePixelOffset:
    """Tests for compute_pixel_offset function."""

    def test_identical_wcs_zero_offset(self):
        """Identical WCS should have zero offset."""
        from astropy.wcs import WCS

        header = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
        }

        wcs = WCS(header)
        dx, dy = compute_pixel_offset(wcs, wcs, 25.0, 25.0)

        assert dx == pytest.approx(0.0, abs=0.01)
        assert dy == pytest.approx(0.0, abs=0.01)

    def test_offset_crpix(self):
        """Different CRPIX should produce pixel offset."""
        from astropy.wcs import WCS

        header_a = {
            'NAXIS': 2,
            'NAXIS1': 100,
            'NAXIS2': 100,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CRPIX1': 50,
            'CRPIX2': 50,
            'CRVAL1': 180.0,
            'CRVAL2': 45.0,
            'CDELT1': -0.001,
            'CDELT2': 0.001,
        }

        header_b = header_a.copy()
        header_b['CRPIX1'] = 55  # Offset by 5 pixels

        wcs_a = WCS(header_a)
        wcs_b = WCS(header_b)

        dx, dy = compute_pixel_offset(wcs_a, wcs_b, 50.0, 50.0)

        # Offset should be approximately 5 pixels in x
        assert abs(dx) == pytest.approx(5.0, abs=0.1)

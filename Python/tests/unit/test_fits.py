"""Tests for heisenberg.io.fits module."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from heisenberg.io.fits import (
    read_fits,
    read_fits_data,
    read_fits_header,
    write_fits,
    get_pixel_scale,
    get_pixel_scale_pc,
    load_tracer_maps,
)

# Skip all tests if astropy is not available
astropy = pytest.importorskip("astropy")
from astropy.io import fits
from astropy.wcs import WCS


def create_test_fits(shape=(100, 100), with_wcs=True, pixel_scale_deg=0.001):
    """Create a test FITS file and return its path."""
    data = np.random.randn(*shape).astype(np.float32)

    header = fits.Header()
    header['TESTKEY'] = 'testvalue'

    if with_wcs:
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CRPIX1'] = shape[1] / 2
        header['CRPIX2'] = shape[0] / 2
        header['CRVAL1'] = 180.0
        header['CRVAL2'] = 45.0
        header['CDELT1'] = -pixel_scale_deg
        header['CDELT2'] = pixel_scale_deg
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'

    hdu = fits.PrimaryHDU(data, header=header)

    with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
        path = Path(f.name)

    hdu.writeto(path, overwrite=True)
    return path, data


class TestReadFits:
    """Tests for read_fits function."""

    def test_reads_data(self):
        """Should read image data correctly."""
        path, expected_data = create_test_fits()
        try:
            data, wcs, header = read_fits(path)
            np.testing.assert_array_almost_equal(data, expected_data)
        finally:
            path.unlink()

    def test_reads_wcs(self):
        """Should read WCS when present."""
        path, _ = create_test_fits(with_wcs=True)
        try:
            data, wcs, header = read_fits(path)
            assert wcs is not None
            assert wcs.naxis == 2
        finally:
            path.unlink()

    def test_no_wcs_returns_none_or_default(self):
        """WCS should be None or default when not explicitly set."""
        path, _ = create_test_fits(with_wcs=False)
        try:
            data, wcs, header = read_fits(path)
            # astropy may create a default WCS, check it's not meaningful
            if wcs is not None:
                # Default WCS has no CTYPE set
                assert wcs.wcs.ctype[0] == '' or wcs.wcs.ctype[0] is None
        finally:
            path.unlink()

    def test_reads_header(self):
        """Should read header keywords."""
        path, _ = create_test_fits()
        try:
            data, wcs, header = read_fits(path)
            assert 'TESTKEY' in header
            assert header['TESTKEY'] == 'testvalue'
        finally:
            path.unlink()


class TestReadFitsData:
    """Tests for read_fits_data function."""

    def test_returns_only_data(self):
        """Should return only the data array."""
        path, expected_data = create_test_fits()
        try:
            data = read_fits_data(path)
            assert isinstance(data, np.ndarray)
            np.testing.assert_array_almost_equal(data, expected_data)
        finally:
            path.unlink()


class TestReadFitsHeader:
    """Tests for read_fits_header function."""

    def test_returns_header_dict(self):
        """Should return header as dictionary."""
        path, _ = create_test_fits()
        try:
            header = read_fits_header(path)
            assert isinstance(header, dict)
            assert 'TESTKEY' in header
        finally:
            path.unlink()


class TestWriteFits:
    """Tests for write_fits function."""

    def test_writes_data(self):
        """Should write data that can be read back."""
        data = np.random.randn(50, 50).astype(np.float32)

        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            path = Path(f.name)

        try:
            write_fits(data, path, overwrite=True)
            read_data = read_fits_data(path)
            np.testing.assert_array_almost_equal(read_data, data)
        finally:
            path.unlink()

    def test_writes_header(self):
        """Should write custom header keywords."""
        data = np.random.randn(50, 50).astype(np.float32)
        header = {'MYKEY': 'myvalue', 'NUMKEY': 42}

        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            path = Path(f.name)

        try:
            write_fits(data, path, header=header, overwrite=True)
            read_header = read_fits_header(path)
            assert read_header['MYKEY'] == 'myvalue'
            assert read_header['NUMKEY'] == 42
        finally:
            path.unlink()

    def test_writes_wcs(self):
        """Should write WCS when provided."""
        data = np.random.randn(50, 50).astype(np.float32)

        # Create a simple WCS
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [25, 25]
        wcs.wcs.cdelt = [-0.001, 0.001]
        wcs.wcs.crval = [180.0, 45.0]
        wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

        with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as f:
            path = Path(f.name)

        try:
            write_fits(data, path, wcs=wcs, overwrite=True)
            _, read_wcs, _ = read_fits(path)
            assert read_wcs is not None
        finally:
            path.unlink()


class TestGetPixelScale:
    """Tests for get_pixel_scale function."""

    def test_returns_scale_from_cdelt(self):
        """Should return scale from CDELT."""
        path, _ = create_test_fits(pixel_scale_deg=0.001)
        try:
            _, wcs, _ = read_fits(path)
            scale = get_pixel_scale(wcs)
            assert np.isclose(scale, 0.001, rtol=1e-5)
        finally:
            path.unlink()

    def test_raises_for_none_wcs(self):
        """Should raise ValueError for None WCS."""
        with pytest.raises(ValueError, match="No WCS"):
            get_pixel_scale(None)


class TestGetPixelScalePc:
    """Tests for get_pixel_scale_pc function."""

    def test_converts_to_parsecs(self):
        """Should convert to parsecs using distance."""
        path, _ = create_test_fits(pixel_scale_deg=0.001)
        try:
            _, wcs, _ = read_fits(path)
            distance_pc = 1e6  # 1 Mpc

            scale_pc = get_pixel_scale_pc(wcs, distance_pc)

            # 0.001 deg at 1 Mpc = 0.001 * pi/180 * 1e6 pc
            expected = 0.001 * np.pi / 180 * 1e6
            assert np.isclose(scale_pc, expected, rtol=1e-3)
        finally:
            path.unlink()


class TestLoadTracerMaps:
    """Tests for load_tracer_maps function."""

    def test_loads_matching_maps(self):
        """Should load star and gas maps with matching shapes."""
        star_path, star_data = create_test_fits(shape=(64, 64))
        gas_path, gas_data = create_test_fits(shape=(64, 64))

        try:
            star, gas, scale, mask = load_tracer_maps(
                star_path, gas_path, distance_pc=1e6
            )
            np.testing.assert_array_almost_equal(star, star_data)
            np.testing.assert_array_almost_equal(gas, gas_data)
            assert scale > 0
            assert mask is None
        finally:
            star_path.unlink()
            gas_path.unlink()

    def test_raises_for_shape_mismatch(self):
        """Should raise for mismatched shapes."""
        star_path, _ = create_test_fits(shape=(64, 64))
        gas_path, _ = create_test_fits(shape=(32, 32))

        try:
            with pytest.raises(ValueError, match="shape mismatch"):
                load_tracer_maps(star_path, gas_path, distance_pc=1e6)
        finally:
            star_path.unlink()
            gas_path.unlink()

    def test_loads_mask(self):
        """Should load mask when provided."""
        star_path, _ = create_test_fits(shape=(64, 64))
        gas_path, _ = create_test_fits(shape=(64, 64))
        mask_path, mask_data = create_test_fits(shape=(64, 64), with_wcs=False)

        try:
            star, gas, scale, mask = load_tracer_maps(
                star_path, gas_path, distance_pc=1e6, mask_path=mask_path
            )
            assert mask is not None
            np.testing.assert_array_almost_equal(mask, mask_data)
        finally:
            star_path.unlink()
            gas_path.unlink()
            mask_path.unlink()

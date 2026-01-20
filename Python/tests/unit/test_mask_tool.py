"""
Unit tests for mask_tool module.

Tests the DS9 region masking workflow including region parsing,
mask creation, and image masking.
"""

import numpy as np
import pytest
from pathlib import Path
import tempfile

from heisenberg.imaging.mask_tool import (
    MaskToolResult,
    mask_tool,
    create_mask_from_ds9,
    region_to_mask,
    regions_to_mask,
)
from heisenberg.imaging.masking import Mask


class TestRegionToMask:
    """Tests for region_to_mask function."""

    def test_circle_region(self):
        """Circle region should create circular mask."""
        from heisenberg.plotting.ds9_regions import DS9Region, RegionType, CoordinateSystem

        region = DS9Region(
            region_type=RegionType.CIRCLE,
            params={'x': 50, 'y': 50, 'radius': 10},
            coordinate_system=CoordinateSystem.IMAGE,
        )
        shape = (100, 100)
        mask = region_to_mask(region, shape, negative=False)

        assert mask.data.shape == shape
        # Center should be unmasked (value close to 1)
        assert mask.data[49, 49] > 0.5
        # Far corner should be masked (value close to 0)
        assert mask.data[0, 0] < 0.5

    def test_circle_region_negative(self):
        """Negative circle region should mask inside."""
        from heisenberg.plotting.ds9_regions import DS9Region, RegionType, CoordinateSystem

        region = DS9Region(
            region_type=RegionType.CIRCLE,
            params={'x': 50, 'y': 50, 'radius': 10},
            coordinate_system=CoordinateSystem.IMAGE,
        )
        shape = (100, 100)
        mask = region_to_mask(region, shape, negative=True)

        assert mask.data.shape == shape
        # Center should be masked (value close to 0)
        assert mask.data[49, 49] < 0.5
        # Far corner should be unmasked (value close to 1)
        assert mask.data[0, 0] > 0.5

    def test_box_region(self):
        """Box region should create rectangular mask."""
        from heisenberg.plotting.ds9_regions import DS9Region, RegionType, CoordinateSystem

        region = DS9Region(
            region_type=RegionType.BOX,
            params={'x': 50, 'y': 50, 'width': 20, 'height': 20, 'angle': 0},
            coordinate_system=CoordinateSystem.IMAGE,
        )
        shape = (100, 100)
        mask = region_to_mask(region, shape, negative=False)

        assert mask.data.shape == shape
        # Center should be unmasked
        assert mask.data[49, 49] > 0.5

    def test_ellipse_region(self):
        """Ellipse region should create elliptical mask."""
        from heisenberg.plotting.ds9_regions import DS9Region, RegionType, CoordinateSystem

        region = DS9Region(
            region_type=RegionType.ELLIPSE,
            params={'x': 50, 'y': 50, 'semi_major': 15, 'semi_minor': 10, 'angle': 0},
            coordinate_system=CoordinateSystem.IMAGE,
        )
        shape = (100, 100)
        mask = region_to_mask(region, shape, negative=False)

        assert mask.data.shape == shape
        # Center should be unmasked
        assert mask.data[49, 49] > 0.5


class TestRegionsToMask:
    """Tests for regions_to_mask function."""

    def test_empty_regions_positive(self):
        """Empty positive regions should return all-masked."""
        shape = (100, 100)
        mask = regions_to_mask([], shape, negative=False)

        assert mask.data.shape == shape
        # All zeros for empty positive mask
        assert np.all(mask.data == 0)

    def test_empty_regions_negative(self):
        """Empty negative regions should return all-unmasked."""
        shape = (100, 100)
        mask = regions_to_mask([], shape, negative=True)

        assert mask.data.shape == shape
        # All ones for empty negative mask
        assert np.all(mask.data == 1)

    def test_multiple_positive_regions(self):
        """Multiple positive regions should be combined with OR."""
        from heisenberg.plotting.ds9_regions import DS9Region, RegionType, CoordinateSystem

        regions = [
            DS9Region(
                region_type=RegionType.CIRCLE,
                params={'x': 25, 'y': 25, 'radius': 10},
                coordinate_system=CoordinateSystem.IMAGE,
            ),
            DS9Region(
                region_type=RegionType.CIRCLE,
                params={'x': 75, 'y': 75, 'radius': 10},
                coordinate_system=CoordinateSystem.IMAGE,
            ),
        ]
        shape = (100, 100)
        mask = regions_to_mask(regions, shape, negative=False)

        # Both circle centers should be unmasked
        assert mask.data[24, 24] > 0.5
        assert mask.data[74, 74] > 0.5
        # Center of image (between circles) should be masked
        assert mask.data[50, 50] < 0.5


class TestMaskTool:
    """Tests for mask_tool function."""

    def test_run_without_masks(self):
        """Should return unchanged image when run_without_masks=True."""
        image = np.random.rand(50, 50)
        result = mask_tool(
            image_input=image,
            run_without_masks=True,
        )

        assert isinstance(result, MaskToolResult)
        np.testing.assert_array_equal(result.masked_image, image)
        assert np.all(result.mask.data == 1)

    def test_raises_without_masks(self):
        """Should raise when no masks provided and run_without_masks=False."""
        image = np.random.rand(50, 50)

        with pytest.raises(ValueError, match="Neither positive nor negative"):
            mask_tool(image_input=image, run_without_masks=False)

    def test_with_ds9_file(self):
        """Should apply mask from DS9 region file."""
        image = np.ones((100, 100)) * 100.0

        # Create a temporary DS9 region file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.reg', delete=False) as f:
            f.write("# Region file format: DS9 version 4.1\n")
            f.write("image\n")
            f.write("circle(50,50,20)\n")
            region_path = Path(f.name)

        try:
            result = mask_tool(
                image_input=image,
                ds9_positive_path=region_path,
            )

            assert isinstance(result, MaskToolResult)
            # Center should be unmasked (kept)
            assert not np.isnan(result.masked_image[49, 49])
            # Corner should be masked (NaN)
            assert np.isnan(result.masked_image[0, 0])
        finally:
            region_path.unlink()

    def test_negative_mask(self):
        """Negative mask should block regions."""
        image = np.ones((100, 100)) * 100.0

        # Create a temporary DS9 region file for negative mask
        with tempfile.NamedTemporaryFile(mode='w', suffix='.reg', delete=False) as f:
            f.write("# Region file format: DS9 version 4.1\n")
            f.write("image\n")
            f.write("circle(50,50,20)\n")
            region_path = Path(f.name)

        try:
            result = mask_tool(
                image_input=image,
                ds9_negative_path=region_path,
            )

            assert isinstance(result, MaskToolResult)
            # Center should be masked (NaN) - blocked by negative mask
            assert np.isnan(result.masked_image[49, 49])
            # Corner should be unmasked (kept)
            assert not np.isnan(result.masked_image[0, 0])
        finally:
            region_path.unlink()


class TestCreateMaskFromDS9:
    """Tests for create_mask_from_ds9 convenience function."""

    def test_creates_mask_from_file(self):
        """Should create mask from DS9 file."""
        shape = (100, 100)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.reg', delete=False) as f:
            f.write("# Region file format: DS9 version 4.1\n")
            f.write("image\n")
            f.write("circle(50,50,15)\n")
            region_path = Path(f.name)

        try:
            mask = create_mask_from_ds9(region_path, shape, negative=False)

            assert isinstance(mask, Mask)
            assert mask.data.shape == shape
            # Center should be unmasked
            assert mask.data[49, 49] > 0.5
        finally:
            region_path.unlink()

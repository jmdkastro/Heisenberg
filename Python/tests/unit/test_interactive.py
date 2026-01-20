"""
Unit tests for interactive peak finding module.

Tests the non-interactive components that can be tested without user input.
"""

import numpy as np
import pytest
from pathlib import Path
import tempfile

from heisenberg.peaks.interactive import (
    InteractivePeakConfig,
    InteractiveResult,
    non_interactive_peak_find,
    _run_peak_detection,
    _print_current_settings,
    _export_region_file,
    _write_report,
)


class TestInteractivePeakConfig:
    """Tests for InteractivePeakConfig dataclass."""

    def test_default_values(self):
        """Should have sensible defaults."""
        config = InteractivePeakConfig()

        assert config.npixmin == 20
        assert config.nsigma == 5.0
        assert config.loglevels is True
        assert config.logrange_s == 2.0
        assert config.logspacing_s == 0.5
        assert config.logrange_g == 2.0
        assert config.logspacing_g == 0.5
        assert config.flux_weighted is False

    def test_custom_values(self):
        """Should accept custom values."""
        config = InteractivePeakConfig(
            npixmin=30,
            nsigma=3.0,
            loglevels=False,
            logrange_s=3.0,
        )

        assert config.npixmin == 30
        assert config.nsigma == 3.0
        assert config.loglevels is False
        assert config.logrange_s == 3.0


class TestInteractiveResult:
    """Tests for InteractiveResult dataclass."""

    def test_create_result(self):
        """Should create result with all fields."""
        config = InteractivePeakConfig()
        result = InteractiveResult(
            star_peaks=[],
            gas_peaks=[],
            config=config,
            finalized=True,
            n_iterations=3,
        )

        assert result.star_peaks == []
        assert result.gas_peaks == []
        assert result.config == config
        assert result.finalized is True
        assert result.n_iterations == 3


class TestRunPeakDetection:
    """Tests for _run_peak_detection helper function."""

    def test_finds_peaks_in_synthetic_image(self):
        """Should find peaks in synthetic Gaussian image."""
        # Create image with clear peaks
        np.random.seed(42)
        image = np.zeros((100, 100))

        # Add Gaussian peaks
        for y, x in [(30, 30), (70, 70), (30, 70)]:
            yy, xx = np.ogrid[:100, :100]
            image += 100 * np.exp(-((xx - x)**2 + (yy - y)**2) / (2 * 5**2))

        # Add noise
        image += np.random.normal(0, 2, image.shape)

        config = InteractivePeakConfig(npixmin=10, nsigma=3.0)
        peaks = _run_peak_detection(image, config, is_star=True)

        # Should find at least the 3 peaks we added
        assert len(peaks) >= 3

    def test_respects_npixmin(self):
        """Higher npixmin should result in fewer peaks."""
        np.random.seed(42)
        image = np.zeros((100, 100))

        # Add several small peaks
        for y, x in [(20, 20), (40, 40), (60, 60), (80, 80)]:
            yy, xx = np.ogrid[:100, :100]
            image += 50 * np.exp(-((xx - x)**2 + (yy - y)**2) / (2 * 3**2))

        image += np.random.normal(0, 1, image.shape)

        config_low = InteractivePeakConfig(npixmin=5, nsigma=2.0)
        config_high = InteractivePeakConfig(npixmin=50, nsigma=2.0)

        peaks_low = _run_peak_detection(image, config_low, is_star=True)
        peaks_high = _run_peak_detection(image, config_high, is_star=True)

        # Higher npixmin should find fewer or equal peaks
        assert len(peaks_high) <= len(peaks_low)


class TestNonInteractivePeakFind:
    """Tests for non_interactive_peak_find function."""

    def test_returns_peaks_for_both_maps(self):
        """Should return star and gas peaks."""
        np.random.seed(42)

        # Create synthetic star and gas images
        star_image = np.zeros((100, 100))
        gas_image = np.zeros((100, 100))

        # Add peaks
        for y, x in [(30, 30), (70, 70)]:
            yy, xx = np.ogrid[:100, :100]
            star_image += 100 * np.exp(-((xx - x)**2 + (yy - y)**2) / (2 * 5**2))
            gas_image += 80 * np.exp(-((xx - x + 5)**2 + (yy - y + 5)**2) / (2 * 5**2))

        star_image += np.random.normal(0, 2, star_image.shape)
        gas_image += np.random.normal(0, 2, gas_image.shape)

        config = InteractivePeakConfig(npixmin=10, nsigma=3.0)
        star_peaks, gas_peaks = non_interactive_peak_find(
            star_image, gas_image, config
        )

        assert len(star_peaks) >= 2
        assert len(gas_peaks) >= 2

    def test_uses_default_config_when_none(self):
        """Should use default config when None provided."""
        np.random.seed(42)
        image = np.random.rand(50, 50) * 10

        star_peaks, gas_peaks = non_interactive_peak_find(image, image, config=None)

        # Should run without error
        assert isinstance(star_peaks, list)
        assert isinstance(gas_peaks, list)


class TestExportRegionFile:
    """Tests for _export_region_file helper function."""

    def test_creates_region_file(self):
        """Should create DS9 region file."""
        # Create mock peaks
        class MockPeak:
            def __init__(self, x, y):
                self.x = x
                self.y = y

        peaks = [MockPeak(10, 20), MockPeak(30, 40)]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "peaks.reg"
            _export_region_file(peaks, output_path, color='red')

            assert output_path.exists()

            # Check file contents
            content = output_path.read_text()
            assert 'point' in content.lower() or 'circle' in content.lower()


class TestWriteReport:
    """Tests for _write_report helper function."""

    def test_creates_report_file(self):
        """Should create report file with parameters."""
        config = InteractivePeakConfig(
            npixmin=25,
            nsigma=4.0,
            logrange_s=2.5,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            report_path = _write_report(config, output_dir, 'NGC300', finalized=True)

            assert report_path.exists()

            content = report_path.read_text()
            assert 'npixmin' in content
            assert '25' in content
            assert 'nsigma' in content
            assert '4.0' in content
            assert 'finalised' in content
            assert '1' in content  # finalized=True -> 1

    def test_report_not_finalized(self):
        """Should indicate not finalized in report."""
        config = InteractivePeakConfig()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            report_path = _write_report(config, output_dir, 'M33', finalized=False)

            content = report_path.read_text()
            assert 'finalised' in content
            # Should have 0 for not finalized
            lines = [l for l in content.split('\n') if 'finalised' in l]
            assert len(lines) == 1
            assert '0' in lines[0]

"""
Tests for configuration loading and parsing.
"""

import pytest
from pathlib import Path
from heisenberg.config import HeisenbergConfig, load_config
from heisenberg.config.reader import parse_input_file, parse_value


class TestParseValue:
    """Tests for value parsing."""

    def test_parse_integer(self):
        """Should parse integers correctly."""
        assert parse_value("42") == 42
        assert parse_value("-17") == -17
        assert parse_value("0") == 0

    def test_parse_float(self):
        """Should parse floats correctly."""
        assert parse_value("3.14") == pytest.approx(3.14)
        assert parse_value("-2.5e-10") == pytest.approx(-2.5e-10)
        assert parse_value("1.0e6") == pytest.approx(1.0e6)
        assert parse_value("840000.") == pytest.approx(840000.0)

    def test_parse_null_marker(self):
        """Should parse '-' as None."""
        assert parse_value("-") is None

    def test_parse_string(self):
        """Should parse strings correctly."""
        assert parse_value("galaxy_name") == "galaxy_name"
        assert parse_value("/path/to/file") == "/path/to/file"


class TestParseInputFile:
    """Tests for input file parsing."""

    def test_parse_sample_file(self, sample_input_file: Path):
        """Should parse the sample input file without errors."""
        params = parse_input_file(sample_input_file)

        # Check some expected keys are present
        assert "datadir" in params
        assert "galaxy" in params
        assert "starfile" in params
        assert "gasfile" in params
        assert "distance" in params

    def test_parse_skips_comments(self, sample_input_file: Path):
        """Should skip comment lines."""
        params = parse_input_file(sample_input_file)

        # Comments should not appear as keys
        assert "#" not in str(list(params.keys()))

    def test_parse_handles_inline_comments(self, sample_input_file: Path):
        """Should strip inline comments from values."""
        params = parse_input_file(sample_input_file)

        # Values should not contain '#' or comment text
        for value in params.values():
            if isinstance(value, str):
                assert "#" not in value


class TestLoadConfig:
    """Tests for full configuration loading."""

    def test_load_sample_config(self, sample_input_file: Path):
        """Should load the sample configuration successfully."""
        config = load_config(sample_input_file)

        assert isinstance(config, HeisenbergConfig)
        assert config.files.galaxy is not None
        assert config.basic_map.distance > 0

    def test_config_file_paths(self, sample_input_file: Path):
        """Should correctly parse file paths."""
        config = load_config(sample_input_file)

        assert isinstance(config.files.datadir, Path)
        assert config.files.starfile is not None
        assert config.files.gasfile is not None

    def test_config_flags(self, sample_input_file: Path):
        """Should correctly parse boolean flags."""
        config = load_config(sample_input_file)

        assert isinstance(config.flags1.mask_images, bool)
        assert isinstance(config.flags1.regrid, bool)
        assert isinstance(config.flags4.tophat, bool)

    def test_config_parameters(self, sample_input_file: Path):
        """Should correctly parse numeric parameters."""
        config = load_config(sample_input_file)

        # Check some expected values from the sample file
        assert config.basic_map.distance == 840000.0
        assert config.basic_map.centrex == 701
        assert config.basic_map.centrey == 701
        assert config.aperture.naperture == 9

    def test_config_optional_values(self, sample_input_file: Path):
        """Should handle optional values (marked with '-')."""
        config = load_config(sample_input_file)

        # These are marked as '-' in the sample file
        assert config.files.starfile2 is None
        assert config.files.gasfile2 is None
        assert config.files.starfile3 is None


class TestConfigValidation:
    """Tests for configuration validation."""

    def test_aperture_validation(self, sample_input_file: Path):
        """Should validate aperture parameter relationships."""
        config = load_config(sample_input_file)

        # peak_res should be <= max_res
        assert config.aperture.peak_res <= config.aperture.max_res

        # Both should be < naperture
        assert config.aperture.peak_res < config.aperture.naperture
        assert config.aperture.max_res < config.aperture.naperture

    def test_config_from_file_method(self, sample_input_file: Path):
        """Should be able to load config using class method."""
        config = HeisenbergConfig.from_file(sample_input_file)

        assert isinstance(config, HeisenbergConfig)


class TestConfigConvenienceProperties:
    """Tests for configuration convenience properties."""

    def test_galaxy_property(self, sample_input_file: Path):
        """Should provide galaxy name through convenience property."""
        config = load_config(sample_input_file)
        assert config.galaxy == config.files.galaxy

    def test_distance_property(self, sample_input_file: Path):
        """Should provide distance through convenience property."""
        config = load_config(sample_input_file)
        assert config.distance == config.basic_map.distance

    def test_datadir_property(self, sample_input_file: Path):
        """Should provide data directory through convenience property."""
        config = load_config(sample_input_file)
        assert config.datadir == config.files.datadir

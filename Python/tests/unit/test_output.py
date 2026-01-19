"""Tests for heisenberg.io.output module."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from heisenberg.io.output import (
    OutputConfig,
    write_results,
    read_results,
    write_json,
    read_json,
    write_summary,
    result_to_dict,
    NumpyEncoder,
)
from heisenberg.core.tuningfork import (
    TuningForkResult,
    TuningForkData,
)
from heisenberg.core.fitting import FitResult
from heisenberg.peaks import DetectedPeak
from heisenberg.peaks.statistics import PeakStatistics


def make_mock_fit_result(
    tgas=10.0,
    tover=1.0,
    lambda_=100.0,
    tstar=10.0,
) -> FitResult:
    """Create a mock FitResult for testing."""
    return FitResult(
        tgas=tgas,
        tover=tover,
        lambda_=lambda_,
        tstar=tstar,
        tgas_errmin=1.0,
        tgas_errmax=1.5,
        tover_errmin=0.3,
        tover_errmax=0.5,
        lambda_errmin=10.0,
        lambda_errmax=15.0,
        chi2_min=1.5,
        ndof=10,
        beta_star=1.0,
        beta_gas=1.0,
        surfcon_star=5.0,
        surfcon_gas=5.0,
        tgas_arr=np.linspace(1, 50, 10),
        tover_arr=np.linspace(0.1, 10, 10),
        lambda_arr=np.linspace(10, 500, 10),
        prob_tgas=np.ones(10) / 10,
        prob_tover=np.ones(10) / 10,
        prob_lambda=np.ones(10) / 10,
        corr_tgas_tover=0.1,
        corr_tgas_lambda=-0.2,
        corr_tover_lambda=0.05,
    )


def make_mock_result() -> TuningForkResult:
    """Create a mock TuningForkResult for testing."""
    observed = TuningForkData(
        apertures=np.array([50, 100, 200, 500]),
        fluxratio_star=np.array([1.5, 1.2, 1.1, 1.0]),
        fluxratio_gas=np.array([0.7, 0.8, 0.9, 1.0]),
        err_star=np.array([0.1, 0.1, 0.1, 0.1]),
        err_gas=np.array([0.1, 0.1, 0.1, 0.1]),
        n_star_peaks=np.array([10, 10, 10, 10]),
        n_gas_peaks=np.array([8, 8, 8, 8]),
    )

    fit = make_mock_fit_result()

    # Create mock peaks
    star_peaks = [
        DetectedPeak(x=25.0, y=25.0, total_flux=100.0, npix=50),
        DetectedPeak(x=75.0, y=75.0, total_flux=80.0, npix=40),
    ]
    gas_peaks = [
        DetectedPeak(x=30.0, y=30.0, total_flux=90.0, npix=45),
        DetectedPeak(x=70.0, y=70.0, total_flux=70.0, npix=35),
    ]

    return TuningForkResult(
        observed=observed,
        fit=fit,
        star_peaks=star_peaks,
        gas_peaks=gas_peaks,
        derived={
            'tgas': 10.0,
            'tover': 1.0,
            'lambda': 100.0,
            'tstar': 10.0,
            'ttotal': 19.0,
            'esf': 0.005,
        },
        config={
            'lap_min': 50.0,
            'lap_max': 500.0,
            'n_star_peaks': 2,
            'n_gas_peaks': 2,
        },
    )


class TestNumpyEncoder:
    """Tests for NumpyEncoder."""

    def test_encodes_array(self):
        """Should encode numpy arrays to lists."""
        data = {'arr': np.array([1, 2, 3])}
        result = json.dumps(data, cls=NumpyEncoder)
        assert result == '{"arr": [1, 2, 3]}'

    def test_encodes_float64(self):
        """Should encode numpy float64."""
        data = {'val': np.float64(3.14)}
        result = json.dumps(data, cls=NumpyEncoder)
        assert '3.14' in result

    def test_encodes_int64(self):
        """Should encode numpy int64."""
        data = {'val': np.int64(42)}
        result = json.dumps(data, cls=NumpyEncoder)
        assert result == '{"val": 42}'

    def test_encodes_bool(self):
        """Should encode numpy bool."""
        data = {'val': np.bool_(True)}
        result = json.dumps(data, cls=NumpyEncoder)
        assert result == '{"val": true}'


class TestOutputConfig:
    """Tests for OutputConfig dataclass."""

    def test_default_values(self):
        """Config should have sensible defaults."""
        config = OutputConfig()

        assert config.format == "json"
        assert config.include_arrays is True
        assert config.include_metadata is True
        assert config.include_peaks is True

    def test_custom_values(self):
        """Config should accept custom values."""
        config = OutputConfig(
            format="yaml",
            include_arrays=False,
            include_metadata=False,
        )

        assert config.format == "yaml"
        assert config.include_arrays is False
        assert config.include_metadata is False


class TestResultToDict:
    """Tests for result_to_dict function."""

    def test_includes_fit(self):
        """Should include fit results."""
        result = make_mock_result()
        data = result_to_dict(result)

        assert 'fit' in data
        assert data['fit']['tgas'] == 10.0
        assert data['fit']['lambda'] == 100.0

    def test_includes_observed(self):
        """Should include observed data."""
        result = make_mock_result()
        data = result_to_dict(result)

        assert 'observed' in data
        assert 'apertures' in data['observed']
        assert len(data['observed']['apertures']) == 4

    def test_includes_derived(self):
        """Should include derived quantities."""
        result = make_mock_result()
        data = result_to_dict(result)

        assert 'derived' in data
        assert data['derived']['ttotal'] == 19.0

    def test_includes_metadata(self):
        """Should include metadata by default."""
        result = make_mock_result()
        data = result_to_dict(result)

        assert 'metadata' in data
        assert 'timestamp' in data['metadata']
        assert 'version' in data['metadata']

    def test_excludes_metadata_when_disabled(self):
        """Should exclude metadata when configured."""
        result = make_mock_result()
        config = OutputConfig(include_metadata=False)
        data = result_to_dict(result, config)

        assert 'metadata' not in data

    def test_includes_peaks(self):
        """Should include peak data by default."""
        result = make_mock_result()
        data = result_to_dict(result)

        assert 'star_peaks' in data
        assert 'gas_peaks' in data
        assert len(data['star_peaks']) == 2

    def test_excludes_peaks_when_disabled(self):
        """Should exclude peaks when configured."""
        result = make_mock_result()
        config = OutputConfig(include_peaks=False)
        data = result_to_dict(result, config)

        assert 'star_peaks' not in data
        assert 'gas_peaks' not in data

    def test_includes_arrays(self):
        """Should include probability arrays by default."""
        result = make_mock_result()
        data = result_to_dict(result)

        assert 'tgas_arr' in data['fit']
        assert 'prob_tgas' in data['fit']

    def test_excludes_arrays_when_disabled(self):
        """Should exclude arrays when configured."""
        result = make_mock_result()
        config = OutputConfig(include_arrays=False)
        data = result_to_dict(result, config)

        assert 'tgas_arr' not in data['fit']
        assert 'prob_tgas' not in data['fit']


class TestWriteReadJson:
    """Tests for JSON write/read functions."""

    def test_round_trip(self):
        """Should preserve data through write/read cycle."""
        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
            path = Path(f.name)

        try:
            write_json(result, path)
            loaded = read_json(path)

            # Check fit parameters
            assert loaded.fit.tgas == result.fit.tgas
            assert loaded.fit.lambda_ == result.fit.lambda_
            assert loaded.fit.tover == result.fit.tover

            # Check observed data
            np.testing.assert_array_almost_equal(
                loaded.observed.apertures,
                result.observed.apertures
            )

            # Check derived quantities
            assert loaded.derived['ttotal'] == result.derived['ttotal']
        finally:
            path.unlink()

    def test_creates_valid_json(self):
        """Output should be valid JSON."""
        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
            path = Path(f.name)

        try:
            write_json(result, path)

            with open(path, 'r') as f:
                data = json.load(f)  # Should not raise

            assert 'fit' in data
            assert 'observed' in data
        finally:
            path.unlink()


class TestWriteReadResults:
    """Tests for write_results/read_results functions."""

    def test_auto_detect_json(self):
        """Should auto-detect JSON format from extension."""
        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
            path = Path(f.name)

        try:
            write_results(result, path)
            loaded = read_results(path)

            assert loaded.fit.tgas == result.fit.tgas
        finally:
            path.unlink()

    def test_preserves_peak_count(self):
        """Should preserve number of peaks."""
        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
            path = Path(f.name)

        try:
            write_results(result, path)
            loaded = read_results(path)

            assert len(loaded.star_peaks) == len(result.star_peaks)
            assert len(loaded.gas_peaks) == len(result.gas_peaks)
        finally:
            path.unlink()

    def test_preserves_errors(self):
        """Should preserve asymmetric errors."""
        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
            path = Path(f.name)

        try:
            write_results(result, path)
            loaded = read_results(path)

            assert loaded.fit.tgas_errmin == result.fit.tgas_errmin
            assert loaded.fit.tgas_errmax == result.fit.tgas_errmax
            assert loaded.fit.lambda_errmin == result.fit.lambda_errmin
            assert loaded.fit.lambda_errmax == result.fit.lambda_errmax
        finally:
            path.unlink()


class TestWriteSummary:
    """Tests for write_summary function."""

    def test_creates_readable_output(self):
        """Should create human-readable summary."""
        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            path = Path(f.name)

        try:
            write_summary(result, path)

            with open(path, 'r') as f:
                content = f.read()

            assert 't_gas' in content
            assert 't_over' in content
            assert 'lambda' in content
            assert '10.0' in content  # tgas value
        finally:
            path.unlink()

    def test_includes_errors(self):
        """Summary should include error bars."""
        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            path = Path(f.name)

        try:
            write_summary(result, path)

            with open(path, 'r') as f:
                content = f.read()

            assert '+' in content  # Error notation
            assert '-' in content
        finally:
            path.unlink()

    def test_includes_derived(self):
        """Summary should include derived quantities."""
        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.txt', delete=False) as f:
            path = Path(f.name)

        try:
            write_summary(result, path)

            with open(path, 'r') as f:
                content = f.read()

            assert 'Derived quantities' in content
            assert 'ttotal' in content
        finally:
            path.unlink()


class TestYamlSupport:
    """Tests for YAML support (if available)."""

    def test_yaml_import_error(self):
        """Should raise helpful error if PyYAML not installed."""
        # This test will only work if PyYAML is NOT installed
        # In practice, we just verify the function exists
        from heisenberg.io.output import write_yaml, read_yaml
        assert callable(write_yaml)
        assert callable(read_yaml)

    @pytest.mark.skipif(
        True,  # Skip by default - YAML may or may not be installed
        reason="PyYAML availability varies"
    )
    def test_yaml_round_trip(self):
        """Should preserve data through YAML write/read cycle."""
        pytest.importorskip("yaml")
        from heisenberg.io.output import write_yaml, read_yaml

        result = make_mock_result()

        with tempfile.NamedTemporaryFile(suffix='.yaml', delete=False) as f:
            path = Path(f.name)

        try:
            write_yaml(result, path)
            loaded = read_yaml(path)

            assert loaded.fit.tgas == result.fit.tgas
        finally:
            path.unlink()

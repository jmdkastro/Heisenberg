"""Tests for heisenberg.plotting module."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

# Skip all tests if matplotlib is not available
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use('Agg')  # Non-interactive backend for testing

from heisenberg.plotting.tuningfork_plot import (
    plot_tuningfork,
    plot_tuningfork_residuals,
    plot_parameter_pdfs,
    plot_summary,
    save_figure,
    check_matplotlib,
)
from heisenberg.core.tuningfork import TuningForkResult, TuningForkData
from heisenberg.core.fitting import FitResult
from heisenberg.peaks import DetectedPeak


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
        tgas_arr=np.linspace(1, 50, 20),
        tover_arr=np.linspace(0.1, 10, 20),
        lambda_arr=np.linspace(10, 500, 20),
        prob_tgas=np.exp(-0.5 * ((np.linspace(1, 50, 20) - 10) / 2)**2),
        prob_tover=np.exp(-0.5 * ((np.linspace(0.1, 10, 20) - 1) / 0.3)**2),
        prob_lambda=np.exp(-0.5 * ((np.linspace(10, 500, 20) - 100) / 20)**2),
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
            'vfb': 5.0,
        },
        config={
            'lap_min': 50.0,
            'lap_max': 500.0,
            'n_star_peaks': 2,
            'n_gas_peaks': 2,
        },
    )


class TestCheckMatplotlib:
    """Tests for matplotlib availability check."""

    def test_does_not_raise_when_available(self):
        """Should not raise when matplotlib is available."""
        check_matplotlib()  # Should not raise


class TestPlotTuningfork:
    """Tests for plot_tuningfork function."""

    def test_creates_figure(self):
        """Should create a matplotlib figure."""
        result = make_mock_result()
        fig, ax = plot_tuningfork(result)

        assert fig is not None
        assert ax is not None
        matplotlib.pyplot.close(fig)

    def test_with_custom_ax(self):
        """Should work with provided axes."""
        result = make_mock_result()
        fig, ax = matplotlib.pyplot.subplots()
        _, returned_ax = plot_tuningfork(result, ax=ax)

        assert returned_ax is ax
        matplotlib.pyplot.close(fig)

    def test_show_model_only(self):
        """Should work with data disabled."""
        result = make_mock_result()
        fig, ax = plot_tuningfork(result, show_data=False)

        assert fig is not None
        matplotlib.pyplot.close(fig)

    def test_show_data_only(self):
        """Should work with model disabled."""
        result = make_mock_result()
        fig, ax = plot_tuningfork(result, show_model=False)

        assert fig is not None
        matplotlib.pyplot.close(fig)

    def test_custom_colors(self):
        """Should accept custom colors."""
        result = make_mock_result()
        fig, ax = plot_tuningfork(result, star_color='red', gas_color='blue')

        assert fig is not None
        matplotlib.pyplot.close(fig)

    def test_linear_scale(self):
        """Should work with linear scales."""
        result = make_mock_result()
        fig, ax = plot_tuningfork(result, logx=False, logy=False)

        assert fig is not None
        matplotlib.pyplot.close(fig)

    def test_with_title(self):
        """Should accept custom title."""
        result = make_mock_result()
        fig, ax = plot_tuningfork(result, title='Test Title')

        assert ax.get_title() == 'Test Title'
        matplotlib.pyplot.close(fig)


class TestPlotTuningforkResiduals:
    """Tests for plot_tuningfork_residuals function."""

    def test_creates_figure(self):
        """Should create a matplotlib figure."""
        result = make_mock_result()
        fig, ax = plot_tuningfork_residuals(result)

        assert fig is not None
        assert ax is not None
        matplotlib.pyplot.close(fig)

    def test_with_custom_ax(self):
        """Should work with provided axes."""
        result = make_mock_result()
        fig, ax = matplotlib.pyplot.subplots()
        _, returned_ax = plot_tuningfork_residuals(result, ax=ax)

        assert returned_ax is ax
        matplotlib.pyplot.close(fig)


class TestPlotParameterPdfs:
    """Tests for plot_parameter_pdfs function."""

    def test_creates_three_axes(self):
        """Should create three axes for tgas, tover, lambda."""
        result = make_mock_result()
        fig, (ax_tgas, ax_tover, ax_lambda) = plot_parameter_pdfs(result)

        assert ax_tgas is not None
        assert ax_tover is not None
        assert ax_lambda is not None
        matplotlib.pyplot.close(fig)

    def test_custom_figsize(self):
        """Should accept custom figure size."""
        result = make_mock_result()
        fig, _ = plot_parameter_pdfs(result, figsize=(15, 5))

        assert fig is not None
        matplotlib.pyplot.close(fig)


class TestPlotSummary:
    """Tests for plot_summary function."""

    def test_creates_figure(self):
        """Should create a multi-panel summary figure."""
        result = make_mock_result()
        fig = plot_summary(result)

        assert fig is not None
        matplotlib.pyplot.close(fig)

    def test_custom_figsize(self):
        """Should accept custom figure size."""
        result = make_mock_result()
        fig = plot_summary(result, figsize=(16, 12))

        assert fig is not None
        matplotlib.pyplot.close(fig)

    def test_with_suptitle(self):
        """Should accept super title."""
        result = make_mock_result()
        fig = plot_summary(result, suptitle='Galaxy Analysis')

        assert fig._suptitle is not None
        matplotlib.pyplot.close(fig)


class TestSaveFigure:
    """Tests for save_figure function."""

    def test_saves_png(self):
        """Should save figure as PNG."""
        result = make_mock_result()
        fig, _ = plot_tuningfork(result)

        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            path = Path(f.name)

        try:
            save_figure(fig, str(path))
            assert path.exists()
            assert path.stat().st_size > 0
        finally:
            matplotlib.pyplot.close(fig)
            path.unlink()

    def test_saves_pdf(self):
        """Should save figure as PDF."""
        result = make_mock_result()
        fig, _ = plot_tuningfork(result)

        with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as f:
            path = Path(f.name)

        try:
            save_figure(fig, str(path))
            assert path.exists()
            assert path.stat().st_size > 0
        finally:
            matplotlib.pyplot.close(fig)
            path.unlink()

    def test_custom_dpi(self):
        """Should accept custom DPI."""
        result = make_mock_result()
        fig, _ = plot_tuningfork(result)

        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            path = Path(f.name)

        try:
            save_figure(fig, str(path), dpi=300)
            assert path.exists()
        finally:
            matplotlib.pyplot.close(fig)
            path.unlink()

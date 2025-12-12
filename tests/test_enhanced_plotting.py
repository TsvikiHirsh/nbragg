"""
Tests for enhanced plotting features (query filtering, error maps, 1D plots).
"""
import pytest
import numpy as np
import tempfile
from pathlib import Path
import matplotlib.pyplot as plt


@pytest.fixture(autouse=True)
def close_figures():
    """Automatically close all matplotlib figures after each test."""
    yield
    plt.close('all')


class TestEnhancedParameterMap:
    """Test enhanced plot_parameter_map with query, errors, and 1D plots."""

    @pytest.fixture
    def grouped_data_2d(self, tmp_path):
        """Create 2x2 grouped data for testing."""
        from nbragg import Data

        for i in range(2):
            for j in range(2):
                channels = np.arange(100, 300)
                # Create varying transmission patterns for different groups
                signal_counts = 1000 - channels * (2 + i + j*0.5) + np.random.randint(-10, 10, len(channels))
                signal_counts = np.maximum(signal_counts, 500)

                signal_file = tmp_path / f"signal_x{i}_y{j}.csv"
                with open(signal_file, 'w') as f:
                    f.write("channel,counts\n")
                    for ch, cnt in zip(channels, signal_counts):
                        f.write(f"{ch},{cnt}\n")

                ob_counts = 1000 + np.random.randint(-5, 5, len(channels))
                ob_file = tmp_path / f"ob_x{i}_y{j}.csv"
                with open(ob_file, 'w') as f:
                    f.write("channel,counts\n")
                    for ch, cnt in zip(channels, ob_counts):
                        f.write(f"{ch},{cnt}\n")

        return Data.from_grouped(
            str(tmp_path / "signal_*.csv"),
            str(tmp_path / "ob_*.csv"),
            L=10, tstep=10e-6, verbosity=0
        )

    @pytest.fixture
    def grouped_data_1d(self, tmp_path):
        """Create 1D grouped data for testing."""
        from nbragg import Data

        for i in range(5):
            channels = np.arange(100, 300)
            signal_counts = 1000 - channels * (2 + i*0.3) + np.random.randint(-10, 10, len(channels))
            signal_counts = np.maximum(signal_counts, 500)

            signal_file = tmp_path / f"signal_{i}.csv"
            with open(signal_file, 'w') as f:
                f.write("channel,counts\n")
                for ch, cnt in zip(channels, signal_counts):
                    f.write(f"{ch},{cnt}\n")

            ob_counts = 1000 + np.random.randint(-5, 5, len(channels))
            ob_file = tmp_path / f"ob_{i}.csv"
            with open(ob_file, 'w') as f:
                f.write("channel,counts\n")
                for ch, cnt in zip(channels, ob_counts):
                    f.write(f"{ch},{cnt}\n")

        return Data.from_grouped(
            str(tmp_path / "signal_*.csv"),
            str(tmp_path / "ob_*.csv"),
            L=10, tstep=10e-6, verbosity=0,
            indices=list(range(5))
        )

    @pytest.fixture
    def fitted_result_2d(self, grouped_data_2d):
        """Create fitted 2D grouped result."""
        from nbragg import TransmissionModel, CrossSection, materials

        xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
        model = TransmissionModel(xs)
        return model.fit(grouped_data_2d, n_jobs=2, progress_bar=False, verbose=False)

    @pytest.fixture
    def fitted_result_1d(self, grouped_data_1d):
        """Create fitted 1D grouped result."""
        from nbragg import TransmissionModel, CrossSection, materials

        xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
        model = TransmissionModel(xs)
        return model.fit(grouped_data_1d, n_jobs=2, progress_bar=False, verbose=False)

    def test_plot_parameter_value(self, fitted_result_2d):
        """Test plotting basic parameter values."""
        # Skip if no successful fits
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        # Should not raise
        ax = fitted_result_2d.plot_parameter_map("norm")
        assert ax is not None

    def test_plot_parameter_error(self, fitted_result_2d):
        """Test plotting parameter errors with _err suffix."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        # Plot error of a parameter
        ax = fitted_result_2d.plot_parameter_map("thickness_err")
        assert ax is not None
        assert "thickness_err" in ax.get_title()

    def test_plot_redchi_map(self, fitted_result_2d):
        """Test plotting redchi statistic."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        ax = fitted_result_2d.plot_parameter_map("redchi")
        assert ax is not None
        assert "redchi" in ax.get_title()

    def test_plot_other_statistics(self, fitted_result_2d):
        """Test plotting other fit statistics."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        # Test various statistics
        for stat in ["chisqr", "aic", "bic"]:
            ax = fitted_result_2d.plot_parameter_map(stat)
            assert ax is not None

    def test_query_filtering(self, fitted_result_2d):
        """Test query filtering of results."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        # Query by redchi - should not raise even if all filtered out
        ax = fitted_result_2d.plot_parameter_map("norm", query="redchi < 10")
        assert ax is not None

    def test_complex_query(self, fitted_result_2d):
        """Test complex query expressions."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        # Complex query with multiple conditions
        ax = fitted_result_2d.plot_parameter_map(
            "thickness",
            query="redchi < 100 and norm > 0.5"
        )
        assert ax is not None

    def test_query_with_parameter_values(self, fitted_result_2d):
        """Test query using parameter values."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        # Query using parameter name directly
        ax = fitted_result_2d.plot_parameter_map(
            "norm",
            query="thickness > 0"
        )
        assert ax is not None

    def test_invalid_query_fallback(self, fitted_result_2d, capsys):
        """Test that invalid query doesn't crash, just warns."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        # Invalid query should print warning but not crash
        ax = fitted_result_2d.plot_parameter_map("norm", query="invalid_param > 0")
        captured = capsys.readouterr()
        assert "Warning" in captured.out or "warning" in captured.out.lower()
        assert ax is not None

    def test_1d_line_plot(self, fitted_result_1d):
        """Test 1D line plot."""
        if len(fitted_result_1d.indices) == 0:
            pytest.skip("No successful fits")

        ax = fitted_result_1d.plot_parameter_map("norm", kind='line')
        assert ax is not None
        # Check that it's a line plot (has lines)
        assert len(ax.get_lines()) > 0

    def test_1d_bar_plot(self, fitted_result_1d):
        """Test 1D bar plot."""
        if len(fitted_result_1d.indices) == 0:
            pytest.skip("No successful fits")

        ax = fitted_result_1d.plot_parameter_map("thickness", kind='bar')
        assert ax is not None
        # Check that it's a bar plot (has patches)
        assert len(ax.patches) > 0

    def test_1d_errorbar_plot(self, fitted_result_1d):
        """Test 1D errorbar plot."""
        if len(fitted_result_1d.indices) == 0:
            pytest.skip("No successful fits")

        ax = fitted_result_1d.plot_parameter_map("norm", kind='errorbar')
        assert ax is not None
        # Errorbar plots have containers
        assert len(ax.containers) > 0 or len(ax.get_lines()) > 0

    def test_1d_invalid_kind_raises(self, fitted_result_1d):
        """Test that invalid kind raises error."""
        if len(fitted_result_1d.indices) == 0:
            pytest.skip("No successful fits")

        with pytest.raises(ValueError, match="Unknown kind"):
            fitted_result_1d.plot_parameter_map("norm", kind='invalid')

    def test_1d_with_query_and_kind(self, fitted_result_1d):
        """Test combining query filtering with 1D plot kind."""
        if len(fitted_result_1d.indices) == 0:
            pytest.skip("No successful fits")

        ax = fitted_result_1d.plot_parameter_map(
            "thickness",
            query="redchi < 100",
            kind='errorbar'
        )
        assert ax is not None

    def test_plot_with_custom_colormap(self, fitted_result_2d):
        """Test custom colormap for 2D plots."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        ax = fitted_result_2d.plot_parameter_map("norm", cmap='plasma')
        assert ax is not None
        # Check colorbar exists
        assert len(ax.figure.axes) > 1  # Plot + colorbar

    def test_plot_with_vmin_vmax(self, fitted_result_2d):
        """Test vmin/vmax for color scaling."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        ax = fitted_result_2d.plot_parameter_map("norm", vmin=0.5, vmax=1.5)
        assert ax is not None

    def test_plot_with_custom_title(self, fitted_result_2d):
        """Test custom title."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        custom_title = "My Custom Title"
        ax = fitted_result_2d.plot_parameter_map("norm", title=custom_title)
        assert ax.get_title() == custom_title

    def test_nonexistent_parameter(self, fitted_result_2d):
        """Test plotting nonexistent parameter gives NaN values."""
        if len(fitted_result_2d.indices) == 0:
            pytest.skip("No successful fits")

        # Should not crash, just plot NaN
        ax = fitted_result_2d.plot_parameter_map("nonexistent_param")
        assert ax is not None

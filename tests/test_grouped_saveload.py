"""
Tests for GroupedFitResult save/load functionality.
"""
import pytest
import numpy as np
import tempfile
import os
from pathlib import Path


class TestGroupedResultSaveLoad:
    """Test saving and loading grouped fit results."""

    @pytest.fixture
    def grouped_data_2d(self, tmp_path):
        """Create 2D grouped test data."""
        from nbragg import Data
        import numpy as np

        # Create simple 2x2 grid with realistic channel values and more data points
        for i in range(2):
            for j in range(2):
                # Generate realistic transmission data
                channels = np.arange(100, 200)
                # Simple decreasing transmission pattern
                signal_counts = 1000 - channels * 2 + np.random.randint(-10, 10, len(channels))
                signal_counts = np.maximum(signal_counts, 500)  # Minimum 500 counts

                # Signal
                signal_file = tmp_path / f"signal_x{i}_y{j}.csv"
                signal_data = "channel,counts\n"
                for ch, cnt in zip(channels, signal_counts):
                    signal_data += f"{ch},{cnt}\n"
                signal_file.write_text(signal_data)

                # Openbeam - mostly flat with some noise
                ob_counts = 1000 + np.random.randint(-5, 5, len(channels))
                ob_file = tmp_path / f"ob_x{i}_y{j}.csv"
                ob_data = "channel,counts\n"
                for ch, cnt in zip(channels, ob_counts):
                    ob_data += f"{ch},{cnt}\n"
                ob_file.write_text(ob_data)

        data = Data.from_grouped(
            str(tmp_path / "signal_*.csv"),
            str(tmp_path / "ob_*.csv"),
            L=10, tstep=10e-6, verbosity=0
        )
        return data

    @pytest.fixture
    def grouped_result(self, grouped_data_2d):
        """Create a grouped fit result."""
        from nbragg import TransmissionModel, CrossSection, materials

        xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
        model = TransmissionModel(xs)
        # Fit with minimal setup
        result = model.fit(grouped_data_2d, n_jobs=2, progress_bar=False, verbose=False)
        return result

    def test_save_and_load_compact(self, grouped_result, tmp_path):
        """Test saving and loading in compact mode."""
        from nbragg import GroupedFitResult

        # Save compact
        save_file = tmp_path / "grouped_compact.json"
        grouped_result.save(str(save_file), compact=True)

        # Check file exists
        assert save_file.exists()

        # Load
        loaded = GroupedFitResult.load(str(save_file))

        # Verify structure
        assert loaded.group_shape == grouped_result.group_shape
        assert len(loaded.indices) == len(grouped_result.indices)
        assert set(loaded.indices) == set(grouped_result.indices)

        # Verify parameters
        for idx in grouped_result.indices:
            orig_result = grouped_result[idx]
            load_result = loaded[idx]

            # Check compact flag
            assert hasattr(load_result, 'compact')
            assert load_result.compact == True

            # Check essential data
            assert hasattr(load_result, 'redchi')
            assert hasattr(load_result, 'params')

            # Check parameter values match
            for param_name in orig_result.params:
                assert param_name in load_result.params
                assert np.isclose(load_result.params[param_name].value,
                                orig_result.params[param_name].value, rtol=1e-10)

    def test_save_and_load_full(self, grouped_result, tmp_path):
        """Test saving and loading full results."""
        from nbragg import GroupedFitResult, TransmissionModel

        # Skip if no successful fits
        if len(grouped_result.indices) == 0:
            pytest.skip("No successful fits to test")

        # Save full
        save_file = tmp_path / "grouped_full.json"
        model_file = tmp_path / "grouped_model.json"
        grouped_result.save(str(save_file), compact=False, model_filename=str(model_file))

        # Check files exist
        assert save_file.exists()
        assert model_file.exists()

        # Load
        loaded = GroupedFitResult.load(str(save_file), model_filename=str(model_file))

        # Verify structure
        assert loaded.group_shape == grouped_result.group_shape
        assert len(loaded.indices) == len(grouped_result.indices)

        # Verify full result data
        for idx in grouped_result.indices:
            orig_result = grouped_result[idx]
            load_result = loaded[idx]

            # Check full result attributes
            assert hasattr(load_result, 'redchi')
            assert hasattr(load_result, 'chisqr')
            assert hasattr(load_result, 'success')
            assert hasattr(load_result, 'nfev')

            # Check statistics match
            if orig_result.redchi is not None:
                assert np.isclose(load_result.redchi, orig_result.redchi, rtol=1e-10)

    def test_compact_saves_smaller_file(self, grouped_result, tmp_path):
        """Verify compact mode creates smaller files."""
        # Skip if no successful fits
        if len(grouped_result.indices) == 0:
            pytest.skip("No successful fits to test")

        compact_file = tmp_path / "grouped_compact.json"
        full_file = tmp_path / "grouped_full.json"

        # Save both
        grouped_result.save(str(compact_file), compact=True)
        grouped_result.save(str(full_file), compact=False, model_filename='')

        # Compact should be smaller
        compact_size = compact_file.stat().st_size
        full_size = full_file.stat().st_size

        assert compact_size < full_size

    def test_loaded_compact_can_plot_map(self, grouped_result, tmp_path):
        """Test that loaded compact results can still create parameter maps."""
        from nbragg import GroupedFitResult

        # Skip if no successful fits
        if len(grouped_result.indices) == 0:
            pytest.skip("No successful fits to test")

        # Save and load compact
        save_file = tmp_path / "grouped_compact.json"
        grouped_result.save(str(save_file), compact=True)
        loaded = GroupedFitResult.load(str(save_file))

        # Should be able to plot parameter map
        ax = loaded.plot_parameter_map('norm')
        assert ax is not None

    def test_single_file_contains_all_results(self, grouped_result, tmp_path):
        """Verify all results are in a single JSON file."""
        import json

        save_file = tmp_path / "grouped.json"
        grouped_result.save(str(save_file), compact=True)

        # Load JSON and verify structure
        with open(save_file, 'r') as f:
            data = json.load(f)

        assert 'version' in data
        assert 'class' in data
        assert data['class'] == 'GroupedFitResult'
        assert 'results' in data
        assert 'indices' in data
        assert 'group_shape' in data

        # All results should be in the 'results' dict
        assert len(data['results']) == len(grouped_result.indices)

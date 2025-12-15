"""
Tests for string indices, folder input, and grouped data addition.
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


class TestStringIndices:
    """Test that all indices are stored as strings."""

    @pytest.fixture
    def grouped_data_2d(self, tmp_path):
        """Create 2x2 grouped data for testing."""
        from nbragg import Data

        # Create signal and openbeam files with 2D grid pattern
        for i in range(2):
            for j in range(2):
                channels = np.arange(100, 200)
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

    def test_indices_are_strings(self, grouped_data_2d):
        """Test that all indices are stored as strings."""
        for idx in grouped_data_2d.indices:
            assert isinstance(idx, str), f"Index {idx} should be a string, got {type(idx)}"

    def test_2d_indices_format(self, grouped_data_2d):
        """Test that 2D indices are formatted as "(x,y)" (no spaces)."""
        # Should have indices like "(0,0)", "(0,1)", "(1,0)", "(1,1)"
        assert "(0,0)" in grouped_data_2d.indices
        assert "(0,1)" in grouped_data_2d.indices
        assert "(1,0)" in grouped_data_2d.indices
        assert "(1,1)" in grouped_data_2d.indices

    def test_access_with_tuple(self, grouped_data_2d):
        """Test accessing groups with tuple indices."""
        # Should be able to access with tuple (0, 0)
        grouped_data_2d.plot(index=(0, 0))
        plt.close()

    def test_access_with_string(self, grouped_data_2d):
        """Test accessing groups with string indices."""
        # Should be able to access with string "(0, 0)"
        grouped_data_2d.plot(index="(0, 0)")
        plt.close()

    def test_access_both_methods_equivalent(self, grouped_data_2d):
        """Test that accessing with tuple or string gives same result."""
        # Get table with tuple
        normalized1 = grouped_data_2d._normalize_index((0, 0))
        table1 = grouped_data_2d.groups[normalized1]

        # Get table with string
        normalized2 = grouped_data_2d._normalize_index("(0, 0)")
        table2 = grouped_data_2d.groups[normalized2]

        # Should be the same table
        assert table1.equals(table2)


class TestFolderInput:
    """Test folder input support for from_grouped."""

    def test_load_from_folder(self, tmp_path):
        """Test loading all CSV files from a folder."""
        from nbragg import Data

        # Create folders
        signal_folder = tmp_path / "signal"
        ob_folder = tmp_path / "openbeam"
        signal_folder.mkdir()
        ob_folder.mkdir()

        # Create 3 CSV files in each folder
        for i in range(3):
            channels = np.arange(100, 200)
            signal_counts = 1000 - channels * 2 + np.random.randint(-10, 10, len(channels))
            signal_counts = np.maximum(signal_counts, 500)

            signal_file = signal_folder / f"pixel_{i}.csv"
            with open(signal_file, 'w') as f:
                f.write("channel,counts\n")
                for ch, cnt in zip(channels, signal_counts):
                    f.write(f"{ch},{cnt}\n")

            ob_counts = 1000 + np.random.randint(-5, 5, len(channels))
            ob_file = ob_folder / f"pixel_{i}.csv"
            with open(ob_file, 'w') as f:
                f.write("channel,counts\n")
                for ch, cnt in zip(channels, ob_counts):
                    f.write(f"{ch},{cnt}\n")

        # Load from folders (not glob patterns)
        data = Data.from_grouped(
            str(signal_folder),
            str(ob_folder),
            L=10, tstep=10e-6, verbosity=0
        )

        assert data.is_grouped
        assert len(data.indices) == 3

    def test_folder_extracts_names(self, tmp_path):
        """Test that named ROI files are extracted correctly."""
        from nbragg import Data

        # Create folders with named ROI files
        signal_folder = tmp_path / "signal"
        ob_folder = tmp_path / "openbeam"
        signal_folder.mkdir()
        ob_folder.mkdir()

        roi_names = ["center", "left", "right"]
        for name in roi_names:
            channels = np.arange(100, 200)
            signal_counts = 1000 - channels * 2 + np.random.randint(-10, 10, len(channels))
            signal_counts = np.maximum(signal_counts, 500)

            signal_file = signal_folder / f"{name}.csv"
            with open(signal_file, 'w') as f:
                f.write("channel,counts\n")
                for ch, cnt in zip(channels, signal_counts):
                    f.write(f"{ch},{cnt}\n")

            ob_counts = 1000 + np.random.randint(-5, 5, len(channels))
            ob_file = ob_folder / f"{name}.csv"
            with open(ob_file, 'w') as f:
                f.write("channel,counts\n")
                for ch, cnt in zip(channels, ob_counts):
                    f.write(f"{ch},{cnt}\n")

        # Load from folders
        data = Data.from_grouped(
            str(signal_folder),
            str(ob_folder),
            L=10, tstep=10e-6, verbosity=0
        )

        # Check that names are extracted (may be in different order)
        assert len(data.indices) == 3
        # Indices should be strings containing the roi names
        assert any("center" in idx for idx in data.indices)
        assert any("left" in idx for idx in data.indices)
        assert any("right" in idx for idx in data.indices)


class TestGroupedDataAddition:
    """Test __add__ method for grouped Data objects."""

    @pytest.fixture
    def two_grouped_datasets(self, tmp_path):
        """Create two identical 2x2 grouped datasets."""
        from nbragg import Data

        datasets = []
        for dataset_num in range(2):
            dataset_dir = tmp_path / f"dataset{dataset_num}"
            dataset_dir.mkdir()

            for i in range(2):
                for j in range(2):
                    channels = np.arange(100, 200)
                    signal_counts = 1000 - channels * 2 + np.random.randint(-10, 10, len(channels))
                    signal_counts = np.maximum(signal_counts, 500)

                    signal_file = dataset_dir / f"signal_x{i}_y{j}.csv"
                    with open(signal_file, 'w') as f:
                        f.write("channel,counts\n")
                        for ch, cnt in zip(channels, signal_counts):
                            f.write(f"{ch},{cnt}\n")

                    ob_counts = 1000 + np.random.randint(-5, 5, len(channels))
                    ob_file = dataset_dir / f"ob_x{i}_y{j}.csv"
                    with open(ob_file, 'w') as f:
                        f.write("channel,counts\n")
                        for ch, cnt in zip(channels, ob_counts):
                            f.write(f"{ch},{cnt}\n")

            data = Data.from_grouped(
                str(dataset_dir / "signal_*.csv"),
                str(dataset_dir / "ob_*.csv"),
                L=10, tstep=10e-6, verbosity=0
            )
            datasets.append(data)

        return datasets

    def test_add_grouped_data(self, two_grouped_datasets):
        """Test adding two grouped Data objects."""
        data1, data2 = two_grouped_datasets

        # Add them together
        combined = data1 + data2

        # Check result is grouped
        assert combined.is_grouped

        # Check indices match
        assert set(combined.indices) == set(data1.indices)
        assert len(combined.indices) == len(data1.indices)

    def test_add_grouped_preserves_shape(self, two_grouped_datasets):
        """Test that addition preserves group shape."""
        data1, data2 = two_grouped_datasets

        combined = data1 + data2

        assert combined.group_shape == data1.group_shape

    def test_add_grouped_has_all_groups(self, two_grouped_datasets):
        """Test that combined data has all groups."""
        data1, data2 = two_grouped_datasets

        combined = data1 + data2

        # Check all groups exist
        for idx in combined.indices:
            assert idx in combined.groups

    def test_add_grouped_mismatched_indices_raises(self, tmp_path):
        """Test that adding grouped data with different indices raises error."""
        from nbragg import Data

        # Create 2x2 grid
        for i in range(2):
            for j in range(2):
                channels = np.arange(100, 200)
                signal_counts = 1000 - channels * 2 + np.random.randint(-10, 10, len(channels))
                signal_counts = np.maximum(signal_counts, 500)

                signal_file = tmp_path / f"signal_x{i}_y{j}.csv"
                with open(signal_file, 'w') as f:
                    f.write("channel,counts\n")
                    for ch, cnt in zip(channels, signal_counts):
                        f.write(f"{ch},{cnt}\n")

                ob_file = tmp_path / f"ob_x{i}_y{j}.csv"
                with open(ob_file, 'w') as f:
                    f.write("channel,counts\n")
                    for ch, cnt in zip(channels, 1000 + np.random.randint(-5, 5, len(channels))):
                        f.write(f"{ch},{cnt}\n")

        data1 = Data.from_grouped(
            str(tmp_path / "signal_*.csv"),
            str(tmp_path / "ob_*.csv"),
            L=10, tstep=10e-6, verbosity=0
        )

        # Create 3x3 grid (different shape)
        tmp_path2 = tmp_path / "other"
        tmp_path2.mkdir()

        for i in range(3):
            for j in range(3):
                channels = np.arange(100, 200)
                signal_counts = 1000 - channels * 2 + np.random.randint(-10, 10, len(channels))
                signal_counts = np.maximum(signal_counts, 500)

                signal_file = tmp_path2 / f"signal_x{i}_y{j}.csv"
                with open(signal_file, 'w') as f:
                    f.write("channel,counts\n")
                    for ch, cnt in zip(channels, signal_counts):
                        f.write(f"{ch},{cnt}\n")

                ob_file = tmp_path2 / f"ob_x{i}_y{j}.csv"
                with open(ob_file, 'w') as f:
                    f.write("channel,counts\n")
                    for ch, cnt in zip(channels, 1000 + np.random.randint(-5, 5, len(channels))):
                        f.write(f"{ch},{cnt}\n")

        data2 = Data.from_grouped(
            str(tmp_path2 / "signal_*.csv"),
            str(tmp_path2 / "ob_*.csv"),
            L=10, tstep=10e-6, verbosity=0
        )

        # Should raise ValueError
        with pytest.raises(ValueError, match="different indices"):
            combined = data1 + data2


class Test1DStringIndices:
    """Test string indices for 1D grouped data."""

    def test_1d_indices_are_strings(self, tmp_path):
        """Test that 1D indices are stored as strings."""
        from nbragg import Data

        for i in range(5):
            channels = np.arange(100, 200)
            signal_counts = 1000 - channels * 2 + np.random.randint(-10, 10, len(channels))
            signal_counts = np.maximum(signal_counts, 500)

            signal_file = tmp_path / f"pixel_{i}.csv"
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

        data = Data.from_grouped(
            str(tmp_path / "pixel_*.csv"),
            str(tmp_path / "ob_*.csv"),
            L=10, tstep=10e-6, verbosity=0
        )

        # All indices should be strings
        for idx in data.indices:
            assert isinstance(idx, str)

        # Should be able to parse back to ints
        for idx in data.indices:
            parsed = data._parse_string_index(idx)
            assert isinstance(parsed, int)

    def test_1d_access_with_int_or_string(self, tmp_path):
        """Test accessing 1D groups with int or string."""
        from nbragg import Data

        for i in range(3):
            channels = np.arange(100, 200)
            signal_counts = 1000 - channels * 2 + np.random.randint(-10, 10, len(channels))
            signal_counts = np.maximum(signal_counts, 500)

            signal_file = tmp_path / f"pixel_{i}.csv"
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

        data = Data.from_grouped(
            str(tmp_path / "pixel_*.csv"),
            str(tmp_path / "ob_*.csv"),
            L=10, tstep=10e-6, verbosity=0
        )

        # Access with int
        data.plot(index=0)
        plt.close()

        # Access with string
        data.plot(index="0")
        plt.close()

        # Both should work without error

"""
Tests for grouped data functionality.

Tests loading, accessing, and analyzing grouped/spatially-resolved data.
"""
import pytest
import pandas as pd
import numpy as np
import tempfile
import os
import shutil
from nbragg import Data


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)


@pytest.fixture
def sample_counts_data():
    """Create sample counts data for testing."""
    tof = np.arange(100, 1000, 10)
    counts = np.random.poisson(1000, len(tof))
    err = np.sqrt(counts)
    return pd.DataFrame({"tof": tof, "counts": counts, "err": err})


@pytest.fixture
def create_test_counts_files():
    """Factory to create test counts files."""
    def _create_files(tmpdir, pattern="2d", n_groups=4):
        """
        Create test counts files.

        Parameters
        ----------
        tmpdir : str
            Directory to create files in
        pattern : str
            "2d" for grid, "1d" for array, "named" for named groups
        n_groups : int
            Number of groups to create
        """
        tof = np.arange(100, 500, 10)

        if pattern == "2d":
            # Create 2x2 grid
            indices = [(0, 0), (0, 1), (1, 0), (1, 1)]
            for x, y in indices:
                # Signal file
                counts = np.random.poisson(1000, len(tof))
                df = pd.DataFrame({
                    'tof': tof,
                    'counts': counts,
                    'err': np.sqrt(counts)
                })
                sig_file = os.path.join(tmpdir, f"signal_x{x}_y{y}.csv")
                df.to_csv(sig_file, index=False, header=False)

                # Openbeam file
                counts_ob = np.random.poisson(1500, len(tof))
                df_ob = pd.DataFrame({
                    'tof': tof,
                    'counts': counts_ob,
                    'err': np.sqrt(counts_ob)
                })
                ob_file = os.path.join(tmpdir, f"openbeam_x{x}_y{y}.csv")
                df_ob.to_csv(ob_file, index=False, header=False)

        elif pattern == "1d":
            # Create 1D array
            for i in range(n_groups):
                # Signal file
                counts = np.random.poisson(1000, len(tof))
                df = pd.DataFrame({
                    'tof': tof,
                    'counts': counts,
                    'err': np.sqrt(counts)
                })
                sig_file = os.path.join(tmpdir, f"signal_pixel_{i}.csv")
                df.to_csv(sig_file, index=False, header=False)

                # Openbeam file
                counts_ob = np.random.poisson(1500, len(tof))
                df_ob = pd.DataFrame({
                    'tof': tof,
                    'counts': counts_ob,
                    'err': np.sqrt(counts_ob)
                })
                ob_file = os.path.join(tmpdir, f"openbeam_pixel_{i}.csv")
                df_ob.to_csv(ob_file, index=False, header=False)

        elif pattern == "named":
            # Create named groups in subdirectories to avoid glob conflicts
            sig_dir = os.path.join(tmpdir, "signal")
            ob_dir = os.path.join(tmpdir, "openbeam")
            os.makedirs(sig_dir, exist_ok=True)
            os.makedirs(ob_dir, exist_ok=True)

            names = ["sample1", "sample2", "sample3", "reference"][:n_groups]
            for name in names:
                # Signal file
                counts = np.random.poisson(1000, len(tof))
                df = pd.DataFrame({
                    'tof': tof,
                    'counts': counts,
                    'err': np.sqrt(counts)
                })
                sig_file = os.path.join(sig_dir, f"{name}.csv")
                df.to_csv(sig_file, index=False, header=False)

                # Openbeam file
                counts_ob = np.random.poisson(1500, len(tof))
                df_ob = pd.DataFrame({
                    'tof': tof,
                    'counts': counts_ob,
                    'err': np.sqrt(counts_ob)
                })
                ob_file = os.path.join(ob_dir, f"{name}.csv")
                df_ob.to_csv(ob_file, index=False, header=False)

        return tmpdir

    return _create_files


class TestFromGrouped2D:
    """Test loading 2D grouped data."""

    def test_load_2d_grid(self, temp_dir, create_test_counts_files):
        """Test loading 2x2 grid with automatic coordinate extraction."""
        create_test_counts_files(temp_dir, pattern="2d")

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        assert data.is_grouped == True
        assert len(data.indices) == 4
        assert data.group_shape == (2, 2)

        # Check indices are strings
        assert all(isinstance(idx, str) for idx in data.indices)
        # Check that they represent 2D tuples like "(0, 0)"
        assert "(0, 0)" in data.indices or (0, 0) in [data._parse_string_index(idx) for idx in data.indices]

        # Check groups dict - can access with tuple or string
        assert len(data.groups) == 4
        assert data._normalize_index((0, 0)) in data.groups
        assert data._normalize_index((1, 1)) in data.groups

    def test_2d_grid_data_content(self, temp_dir, create_test_counts_files):
        """Test that loaded 2D data has correct structure."""
        create_test_counts_files(temp_dir, pattern="2d")

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Check each group has proper data
        for idx in data.indices:
            group_table = data.groups[idx]
            assert 'wavelength' in group_table.columns
            assert 'trans' in group_table.columns
            assert 'err' in group_table.columns
            assert len(group_table) > 0

    def test_2d_grid_default_table(self, temp_dir, create_test_counts_files):
        """Test that default table is set to first group."""
        create_test_counts_files(temp_dir, pattern="2d")

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Default table should be first index
        first_idx = data.indices[0]
        assert data.table is data.groups[first_idx]


class TestFromGrouped1D:
    """Test loading 1D grouped data."""

    def test_load_1d_array(self, temp_dir, create_test_counts_files):
        """Test loading 1D array with automatic index extraction."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=5)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        assert data.is_grouped == True
        assert len(data.indices) == 5
        assert data.group_shape == (5,)

        # Check indices are strings (converted from ints)
        assert all(isinstance(idx, str) for idx in data.indices)

    def test_1d_with_custom_indices(self, temp_dir, create_test_counts_files):
        """Test 1D loading with user-provided indices."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=3)

        custom_indices = [10, 20, 30]

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6,
            indices=custom_indices
        )

        # Indices should be converted to strings
        assert data.indices == ["10", "20", "30"]
        assert "10" in data.groups or data._normalize_index(10) in data.groups
        assert "20" in data.groups or data._normalize_index(20) in data.groups
        assert "30" in data.groups or data._normalize_index(30) in data.groups


class TestFromGroupedNamed:
    """Test loading named groups."""

    def test_load_named_groups(self, temp_dir, create_test_counts_files):
        """Test loading with named indices."""
        create_test_counts_files(temp_dir, pattern="named", n_groups=3)

        custom_indices = ["sample1", "sample2", "sample3"]

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal", "*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam", "*.csv"),
            L=10.0,
            tstep=10e-6,
            indices=custom_indices
        )

        assert data.is_grouped == True
        assert len(data.indices) == 3
        assert data.group_shape is None  # Named groups have no shape

        # Check indices are strings
        assert all(isinstance(idx, str) for idx in data.indices)
        assert "sample1" in data.groups


class TestGroupedWithCorrections:
    """Test grouped data with L0, t0, dropna, empty region corrections."""

    def test_grouped_with_l0_t0(self, temp_dir, create_test_counts_files):
        """Test that L0 and t0 corrections apply to all groups."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=2)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6,
            L0=1.05,
            t0=5.0
        )

        # Check that corrections were applied (data should load without error)
        assert len(data.groups) == 2
        for idx in data.indices:
            assert 'wavelength' in data.groups[idx].columns

    def test_grouped_with_dropna(self, temp_dir, create_test_counts_files):
        """Test dropna parameter with grouped data."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=2)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6,
            dropna=True
        )

        # Check no NaN values in any group
        for idx in data.indices:
            assert not data.groups[idx].isna().any().any()


class TestPatternExtraction:
    """Test coordinate extraction from filenames."""

    def test_extract_2d_coordinates(self, temp_dir, create_test_counts_files):
        """Test extraction of x, y coordinates from filenames."""
        create_test_counts_files(temp_dir, pattern="2d")

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Check extracted coordinates (now string format without spaces)
        assert "(0,0)" in data.indices
        assert "(0,1)" in data.indices
        assert "(1,0)" in data.indices
        assert "(1,1)" in data.indices

    def test_extract_1d_indices(self, temp_dir, create_test_counts_files):
        """Test extraction of 1D indices from filenames."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=4)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Check extracted indices (now string format)
        assert "0" in data.indices
        assert "1" in data.indices
        assert "2" in data.indices
        assert "3" in data.indices


class TestGroupedDataAccess:
    """Test accessing grouped data."""

    def test_access_specific_group(self, temp_dir, create_test_counts_files):
        """Test accessing a specific group by index."""
        create_test_counts_files(temp_dir, pattern="2d")

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Access specific group (with normalized index)
        group_0_0 = data.groups[data._normalize_index((0, 0))]
        assert isinstance(group_0_0, pd.DataFrame)
        assert 'wavelength' in group_0_0.columns

    def test_iterate_groups(self, temp_dir, create_test_counts_files):
        """Test iterating over all groups."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=3)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Iterate over groups
        count = 0
        for idx in data.indices:
            group = data.groups[idx]
            assert isinstance(group, pd.DataFrame)
            count += 1

        assert count == 3


class TestErrorHandling:
    """Test error handling for grouped data."""

    def test_no_files_found(self, temp_dir):
        """Test error when no files match pattern."""
        with pytest.raises(ValueError, match="No files found"):
            Data.from_grouped(
                signal=os.path.join(temp_dir, "nonexistent_*.csv"),
                openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
                L=10.0,
                tstep=10e-6
            )

    def test_mismatch_file_count(self, temp_dir, create_test_counts_files):
        """Test error when signal and openbeam file counts don't match."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=3)

        # Remove one openbeam file
        ob_files = [f for f in os.listdir(temp_dir) if 'openbeam' in f]
        os.remove(os.path.join(temp_dir, ob_files[0]))

        with pytest.raises(ValueError, match="Mismatch"):
            Data.from_grouped(
                signal=os.path.join(temp_dir, "signal_*.csv"),
                openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
                L=10.0,
                tstep=10e-6
            )

    def test_mismatch_indices_count(self, temp_dir, create_test_counts_files):
        """Test error when custom indices count doesn't match files."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=3)

        with pytest.raises(ValueError, match="Number of indices"):
            Data.from_grouped(
                signal=os.path.join(temp_dir, "signal_*.csv"),
                openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
                L=10.0,
                tstep=10e-6,
                indices=[0, 1]  # Only 2 indices for 3 files
            )


class TestGroupedPlotting:
    """Test plotting functionality for grouped data."""

    def test_plot_specific_group_2d(self, temp_dir, create_test_counts_files):
        """Test plotting a specific group from 2D data."""
        create_test_counts_files(temp_dir, pattern="2d")

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Plot specific group - should not raise
        ax = data.plot(index=(0, 1))
        assert ax is not None

    def test_plot_specific_group_1d(self, temp_dir, create_test_counts_files):
        """Test plotting a specific group from 1D data."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=3)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Plot specific group - should not raise
        ax = data.plot(index=1)
        assert ax is not None

    def test_plot_default_group(self, temp_dir, create_test_counts_files):
        """Test plotting default (first) group when no index specified."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=3)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Plot without index - should plot first group
        ax = data.plot()
        assert ax is not None

    def test_plot_invalid_index(self, temp_dir, create_test_counts_files):
        """Test error when plotting with invalid index."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=3)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Try to plot non-existent index
        with pytest.raises(ValueError, match="Index.*not found"):
            data.plot(index=99)

    def test_plot_index_on_non_grouped(self, sample_counts_data):
        """Test error when trying to use index on non-grouped data."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1

        data = Data.from_counts(signal, openbeam, L=10, tstep=10e-6)

        # Should raise error when using index on non-grouped data
        with pytest.raises(ValueError, match="Cannot specify index for non-grouped data"):
            data.plot(index=0)


class TestPlotMap:
    """Test plot_map functionality for visualizing average transmission."""

    def test_plot_map_2d(self, temp_dir, create_test_counts_files):
        """Test plot_map for 2D grid data."""
        create_test_counts_files(temp_dir, pattern="2d")

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Plot map - should not raise
        ax = data.plot_map(wlmin=1.0, wlmax=3.0)
        assert ax is not None

    def test_plot_map_1d(self, temp_dir, create_test_counts_files):
        """Test plot_map for 1D array data."""
        create_test_counts_files(temp_dir, pattern="1d", n_groups=4)

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Plot map - should not raise
        ax = data.plot_map(wlmin=1.0, wlmax=3.0)
        assert ax is not None

    def test_plot_map_named(self, temp_dir, create_test_counts_files):
        """Test plot_map for named groups."""
        create_test_counts_files(temp_dir, pattern="named", n_groups=3)

        custom_indices = ["sample1", "sample2", "sample3"]

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal", "*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam", "*.csv"),
            L=10.0,
            tstep=10e-6,
            indices=custom_indices
        )

        # Plot map - should create bar chart
        ax = data.plot_map(wlmin=1.0, wlmax=3.0)
        assert ax is not None

    def test_plot_map_with_custom_params(self, temp_dir, create_test_counts_files):
        """Test plot_map with custom parameters."""
        create_test_counts_files(temp_dir, pattern="2d")

        data = Data.from_grouped(
            signal=os.path.join(temp_dir, "signal_*.csv"),
            openbeam=os.path.join(temp_dir, "openbeam_*.csv"),
            L=10.0,
            tstep=10e-6
        )

        # Plot with custom parameters
        ax = data.plot_map(wlmin=1.5, wlmax=4.0, cmap='plasma', title='Custom Title')
        assert ax is not None

    def test_plot_map_on_non_grouped_raises(self, sample_counts_data):
        """Test that plot_map raises error on non-grouped data."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1

        data = Data.from_counts(signal, openbeam, L=10, tstep=10e-6)

        # Should raise error
        with pytest.raises(ValueError, match="plot_map only works for grouped data"):
            data.plot_map()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

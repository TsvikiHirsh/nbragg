"""
Tests for Data class improvements including:
- L0 and t0 corrections
- DataFrame input support
- dropna functionality
- sys_err removal
"""
import pytest
import pandas as pd
import numpy as np
import nbragg
from nbragg import Data
import tempfile
import os


@pytest.fixture
def sample_counts_data():
    """Create sample counts data for testing."""
    tof = np.arange(100, 1000, 10)
    counts = np.random.poisson(1000, len(tof))
    err = np.sqrt(counts)
    return pd.DataFrame({"tof": tof, "counts": counts, "err": err})


@pytest.fixture
def sample_transmission_data():
    """Create sample transmission data for testing."""
    wavelength = np.linspace(0.5, 5.0, 50)
    trans = np.exp(-0.1 * wavelength)
    err = 0.01 * trans
    return pd.DataFrame({"wavelength": wavelength, "trans": trans, "err": err})


class TestDataFrameInput:
    """Test that Data methods accept pandas DataFrames."""

    def test_from_counts_with_dataframes(self, sample_counts_data):
        """Test from_counts accepts DataFrames for signal and openbeam."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1  # Make openbeam slightly different

        data = Data.from_counts(signal, openbeam, L=10, tstep=10e-6)

        assert data.table is not None
        assert "wavelength" in data.table.columns
        assert "trans" in data.table.columns
        assert "err" in data.table.columns
        assert len(data.table) == len(signal)

    def test_from_counts_with_empty_dataframes(self, sample_counts_data):
        """Test from_counts accepts DataFrames for empty region correction."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1

        empty_signal = sample_counts_data.copy()
        empty_signal["counts"] *= 0.1
        empty_openbeam = sample_counts_data.copy()
        empty_openbeam["counts"] *= 0.1

        data = Data.from_counts(
            signal, openbeam,
            empty_signal=empty_signal,
            empty_openbeam=empty_openbeam,
            L=10, tstep=10e-6
        )

        assert data.table is not None
        assert len(data.table) == len(signal)

    def test_from_transmission_with_dataframe(self, sample_transmission_data):
        """Test from_transmission accepts DataFrames."""
        data = Data.from_transmission(sample_transmission_data)

        assert data.table is not None
        assert "wavelength" in data.table.columns
        assert "trans" in data.table.columns
        assert "err" in data.table.columns
        assert len(data.table) == len(sample_transmission_data)


class TestL0T0Corrections:
    """Test L0 and t0 correction parameters."""

    def test_l0_correction_shifts_wavelength(self, sample_counts_data):
        """Test that L0 correction shifts wavelength values.

        L0 is a scale factor (default 1.0):
        - L0 > 1.0: longer flight path, TOF is reduced (dtof negative)
        - L0 < 1.0: shorter flight path, TOF is increased (dtof positive)
        """
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1

        L = 10.0
        L0 = 1.1  # 10% longer path
        tstep = 10e-6

        # Create data without correction (L0=1.0, default)
        data_no_correction = Data.from_counts(signal, openbeam, L=L, tstep=tstep)

        # Create data with L0 correction
        data_with_correction = Data.from_counts(signal, openbeam, L=L, L0=L0, tstep=tstep)

        # With L0 = 1.1, dtof = (1.0 - 1.1)*tof + 0 = -0.1*tof
        # corrected_tof = tof + dtof = 0.9*tof (shorter)
        # Shorter TOF → higher energy → shorter wavelength
        assert not np.allclose(
            data_no_correction.table["wavelength"].values,
            data_with_correction.table["wavelength"].values
        )
        # Check that wavelengths are indeed smaller with L0 > 1.0
        assert np.all(
            data_with_correction.table["wavelength"].values <
            data_no_correction.table["wavelength"].values
        )

    def test_t0_correction_shifts_wavelength(self, sample_counts_data):
        """Test that t0 correction shifts wavelength values."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1

        L = 10.0
        t0 = 10  # Positive offset in tof units
        tstep = 10e-6

        # Create data without correction
        data_no_correction = Data.from_counts(signal, openbeam, L=L, tstep=tstep)

        # Create data with t0 correction
        data_with_correction = Data.from_counts(signal, openbeam, L=L, t0=t0, tstep=tstep)

        # With t0 > 0, effective time is larger, so wavelength should be larger
        assert not np.allclose(
            data_no_correction.table["wavelength"].values,
            data_with_correction.table["wavelength"].values
        )
        # Check that wavelengths are indeed larger with positive t0
        assert np.all(
            data_with_correction.table["wavelength"].values >
            data_no_correction.table["wavelength"].values
        )

    def test_combined_l0_t0_corrections(self, sample_counts_data):
        """Test that L0 and t0 can be applied together."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1

        data = Data.from_counts(
            signal, openbeam,
            L=10.0, L0=0.5, t0=10,
            tstep=10e-6
        )

        assert data.table is not None
        assert len(data.table) == len(signal)

    def test_l0_less_than_one(self, sample_counts_data):
        """Test that L0 < 1.0 works (shorter flight path)."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1

        # L0 < 1.0 means shorter flight path
        data = Data.from_counts(
            signal, openbeam,
            L=10.0, L0=0.9, t0=-5,
            tstep=10e-6
        )

        assert data.table is not None
        assert len(data.table) == len(signal)


class TestDropNA:
    """Test dropna functionality."""

    def test_dropna_parameter_in_from_counts(self, sample_counts_data):
        """Test dropna parameter in from_counts."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()
        openbeam["counts"] *= 1.1

        # Introduce some NaN values by dividing by zero
        signal.loc[5, "counts"] = 0
        openbeam.loc[5, "counts"] = 0

        # Without dropna
        data_with_nan = Data.from_counts(signal, openbeam, L=10, tstep=10e-6, dropna=False)

        # With dropna
        data_without_nan = Data.from_counts(signal, openbeam, L=10, tstep=10e-6, dropna=True)

        # Data with dropna should have fewer rows
        assert len(data_without_nan.table) < len(data_with_nan.table)
        # Data without dropna should have NaN values
        assert data_with_nan.table.isna().any().any()
        # Data with dropna should have no NaN values
        assert not data_without_nan.table.isna().any().any()

    def test_dropna_parameter_in_from_transmission(self, sample_transmission_data):
        """Test dropna parameter in from_transmission."""
        df = sample_transmission_data.copy()
        # Introduce NaN
        df.loc[5, "trans"] = np.nan

        # Without dropna
        data_with_nan = Data.from_transmission(df, dropna=False)

        # With dropna
        data_without_nan = Data.from_transmission(df, dropna=True)

        # Data with dropna should have fewer rows
        assert len(data_without_nan.table) < len(data_with_nan.table)
        assert not data_without_nan.table.isna().any().any()

    def test_dropna_method_inplace_false(self, sample_transmission_data):
        """Test dropna method with inplace=False."""
        df = sample_transmission_data.copy()
        df.loc[5, "trans"] = np.nan

        data = Data.from_transmission(df, dropna=False)

        # Call dropna with inplace=False
        data_clean = data.dropna(inplace=False)

        # Original should still have NaN
        assert data.table.isna().any().any()
        # New object should not have NaN
        assert not data_clean.table.isna().any().any()
        # Should be different objects
        assert data is not data_clean

    def test_dropna_method_inplace_true(self, sample_transmission_data):
        """Test dropna method with inplace=True."""
        df = sample_transmission_data.copy()
        df.loc[5, "trans"] = np.nan

        data = Data.from_transmission(df, dropna=False)
        original_id = id(data)

        # Call dropna with inplace=True
        result = data.dropna(inplace=True)

        # Should return same object
        assert result is data
        assert id(result) == original_id
        # Should not have NaN anymore
        assert not data.table.isna().any().any()

    def test_dropna_method_usage_pattern(self, sample_transmission_data):
        """Test the usage pattern data = data.dropna()."""
        df = sample_transmission_data.copy()
        df.loc[5, "trans"] = np.nan

        data = Data.from_transmission(df, dropna=False)

        # Usage pattern: data = data.dropna()
        data = data.dropna()

        assert not data.table.isna().any().any()


class TestSysErrRemoval:
    """Test that sys_err has been properly removed."""

    def test_from_counts_no_sys_err_parameter(self, sample_counts_data):
        """Test that from_counts doesn't accept sys_err parameter."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()

        # This should raise TypeError because sys_err is no longer a parameter
        with pytest.raises(TypeError, match="sys_err"):
            Data.from_counts(signal, openbeam, L=10, tstep=10e-6, sys_err=0.01)

    def test_data_object_no_sys_err_attribute(self, sample_counts_data):
        """Test that Data objects don't have sys_err attribute after initialization."""
        signal = sample_counts_data.copy()
        openbeam = sample_counts_data.copy()

        data = Data.from_counts(signal, openbeam, L=10, tstep=10e-6)

        # sys_err should not be an attribute
        assert not hasattr(data, 'sys_err') or data.sys_err is None


class TestUncertaintyCalculations:
    """Test that uncertainty calculations are correct."""

    def test_basic_transmission_uncertainty(self):
        """Test basic transmission uncertainty formula."""
        # Create simple data where we can verify the calculation
        signal_counts = np.array([1000.0] * 10)
        openbeam_counts = np.array([2000.0] * 10)
        signal_err = np.sqrt(signal_counts)
        openbeam_err = np.sqrt(openbeam_counts)

        signal = pd.DataFrame({
            "tof": np.arange(100, 200, 10),
            "counts": signal_counts,
            "err": signal_err
        })

        openbeam = pd.DataFrame({
            "tof": np.arange(100, 200, 10),
            "counts": openbeam_counts,
            "err": openbeam_err
        })

        data = Data.from_counts(signal, openbeam, L=10, tstep=10e-6)

        # Expected transmission
        expected_trans = signal_counts / openbeam_counts

        # Expected error: trans * sqrt((σ_S/S)² + (σ_O/O)²)
        expected_err = expected_trans * np.sqrt(
            (signal_err / signal_counts)**2 +
            (openbeam_err / openbeam_counts)**2
        )

        np.testing.assert_allclose(
            data.table["trans"].values,
            expected_trans,
            rtol=1e-10
        )

        np.testing.assert_allclose(
            data.table["err"].values,
            expected_err,
            rtol=1e-10
        )

    def test_empty_region_corrected_uncertainty(self):
        """Test uncertainty calculation with empty region correction."""
        # Create simple data
        signal_counts = np.array([1000.0] * 10)
        openbeam_counts = np.array([2000.0] * 10)
        empty_signal_counts = np.array([100.0] * 10)
        empty_openbeam_counts = np.array([200.0] * 10)

        signal = pd.DataFrame({
            "tof": np.arange(100, 200, 10),
            "counts": signal_counts,
            "err": np.sqrt(signal_counts)
        })

        openbeam = pd.DataFrame({
            "tof": np.arange(100, 200, 10),
            "counts": openbeam_counts,
            "err": np.sqrt(openbeam_counts)
        })

        empty_signal = pd.DataFrame({
            "tof": np.arange(100, 200, 10),
            "counts": empty_signal_counts,
            "err": np.sqrt(empty_signal_counts)
        })

        empty_openbeam = pd.DataFrame({
            "tof": np.arange(100, 200, 10),
            "counts": empty_openbeam_counts,
            "err": np.sqrt(empty_openbeam_counts)
        })

        data = Data.from_counts(
            signal, openbeam,
            empty_signal=empty_signal,
            empty_openbeam=empty_openbeam,
            L=10, tstep=10e-6
        )

        # Expected transmission
        expected_trans = (signal_counts / openbeam_counts) * (empty_openbeam_counts / empty_signal_counts)

        # Expected error: trans * sqrt((σ_S/S)² + (σ_O/O)² + (σ_ES/ES)² + (σ_EO/EO)²)
        expected_err = expected_trans * np.sqrt(
            (np.sqrt(signal_counts) / signal_counts)**2 +
            (np.sqrt(openbeam_counts) / openbeam_counts)**2 +
            (np.sqrt(empty_signal_counts) / empty_signal_counts)**2 +
            (np.sqrt(empty_openbeam_counts) / empty_openbeam_counts)**2
        )

        np.testing.assert_allclose(
            data.table["trans"].values,
            expected_trans,
            rtol=1e-10
        )

        np.testing.assert_allclose(
            data.table["err"].values,
            expected_err,
            rtol=1e-10
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

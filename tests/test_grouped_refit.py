"""
Tests for GroupedFitResult.fit() — selective group refitting.
"""
import pytest
import numpy as np
import pandas as pd
import warnings

from nbragg import CrossSection, materials
from nbragg.models import TransmissionModel
from nbragg.data import Data
from nbragg.grouped_fit import GroupedFitResult


def _make_grouped_data(model, n_groups=4, n_points=80, noise_scale=0.005):
    """
    Create synthetic grouped Data with known transmission curves.

    Returns Data with 1D group_shape=(n_groups,).
    """
    wl = np.linspace(1.5, 5.0, n_points)
    data = Data()
    data.is_grouped = True
    data.groups = {}
    data.indices = []
    data.group_shape = (n_groups,)
    data.L = None
    data.tstep = None

    for i in range(n_groups):
        # Slightly different thickness per group
        thickness = 0.5 + 0.1 * i
        trans = model.transmission(wl, thickness=thickness, norm=1.0)
        noise = np.random.default_rng(seed=42 + i).normal(0, noise_scale, n_points)
        table = pd.DataFrame({
            'wavelength': wl,
            'trans': trans + noise,
            'err': np.full(n_points, noise_scale),
        })
        idx = str(i)
        data.groups[idx] = table
        data.indices.append(idx)

    return data


@pytest.fixture(scope="module")
def model():
    """Shared model fixture."""
    xs = CrossSection(iron=materials["Fe_sg229_Iron-alpha"])
    return TransmissionModel(xs)


@pytest.fixture(scope="module")
def grouped_data(model):
    """Shared synthetic grouped data fixture."""
    return _make_grouped_data(model, n_groups=4)


@pytest.fixture(scope="module")
def initial_result(model, grouped_data):
    """
    Fit all groups once (sequential, fast) so refit tests can work from this.
    """
    result = model.fit(
        grouped_data,
        backend="sequential",
        wlmin=1.5,
        wlmax=5.0,
        method="leastsq",
        progress_bar=False,
        verbose=False,
    )
    return result


class TestGroupedRefit:
    """Tests for GroupedFitResult.fit()."""

    def test_initial_result_has_model(self, initial_result):
        """Model reference is stored after grouped fit."""
        assert initial_result.model is not None
        assert isinstance(initial_result.model, TransmissionModel)

    def test_refit_with_query_sequential(self, initial_result, grouped_data):
        """Refit groups matching a query via sequential backend."""
        # Use a query that selects some groups (redchi > 0 selects all, but we
        # just need to verify the machinery works)
        result2 = initial_result.fit(
            grouped_data,
            query="redchi > 0",
            backend="sequential",
            wlmin=1.5,
            wlmax=5.0,
            method="leastsq",
            progress_bar=False,
        )
        # All groups present in merged result
        assert len(result2) == len(initial_result)
        assert set(result2.indices) == set(initial_result.indices)
        assert result2.group_shape == initial_result.group_shape
        assert result2.model is not None

    def test_refit_with_explicit_indices(self, initial_result, grouped_data):
        """Refit specific groups by explicit index list."""
        # Refit only group 0 and group 2
        result2 = initial_result.fit(
            grouped_data,
            indices=[0, 2],
            backend="sequential",
            wlmin=1.5,
            wlmax=5.0,
            method="leastsq",
            progress_bar=False,
        )
        # All groups present
        assert len(result2) == len(initial_result)
        # Non-refitted groups should be identical objects (copied)
        assert result2["1"] is initial_result["1"]
        assert result2["3"] is initial_result["3"]
        # Refitted groups should be different objects
        assert result2["0"] is not initial_result["0"]
        assert result2["2"] is not initial_result["2"]

    def test_refit_uses_previous_params_as_init(self, initial_result, grouped_data):
        """Refitted groups use previous best-fit params as starting values."""
        # Get previous best-fit thickness for group 0
        prev_thickness = initial_result["0"].params["thickness"].value

        # Refit group 0
        result2 = initial_result.fit(
            grouped_data,
            indices=[0],
            backend="sequential",
            wlmin=1.5,
            wlmax=5.0,
            method="leastsq",
            progress_bar=False,
        )

        # The refitted thickness should be very close to previous
        # (same data, same starting point → converges to same answer)
        new_thickness = result2["0"].params["thickness"].value
        assert abs(new_thickness - prev_thickness) < 0.01

    def test_refit_with_params_override(self, initial_result, grouped_data):
        """params argument overrides starting values for refitted groups."""
        override_params = initial_result.model.params.copy()
        override_params["thickness"].value = 0.3  # Different starting point

        result2 = initial_result.fit(
            grouped_data,
            indices=[0],
            params=override_params,
            backend="sequential",
            wlmin=1.5,
            wlmax=5.0,
            method="leastsq",
            progress_bar=False,
        )

        # Should still converge to a reasonable answer (the fit is well-constrained)
        assert result2["0"].params["thickness"].value > 0

    def test_empty_query_returns_copy_with_warning(self, initial_result, grouped_data):
        """Query matching zero groups warns and returns copy of original."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result2 = initial_result.fit(
                grouped_data,
                query="redchi < -999",
                backend="sequential",
                progress_bar=False,
            )
            # Should have raised a warning
            assert any("zero groups" in str(warning.message).lower() for warning in w)

        # All groups should be present and identical to original
        assert len(result2) == len(initial_result)
        for idx in result2.indices:
            assert result2[idx] is initial_result[idx]

    def test_no_model_raises_error(self):
        """fit() on result with model=None raises ValueError."""
        result = GroupedFitResult(group_shape=(2,))
        assert result.model is None

        with pytest.raises(ValueError, match="No model reference"):
            result.fit(None, query="redchi > 0")

    def test_no_query_or_indices_raises_error(self, initial_result, grouped_data):
        """Must provide at least one of query or indices."""
        with pytest.raises(ValueError, match="Must provide at least one"):
            initial_result.fit(grouped_data)

    def test_invalid_index_raises_error(self, initial_result, grouped_data):
        """Invalid index in indices list raises KeyError."""
        with pytest.raises(KeyError, match="not found"):
            initial_result.fit(
                grouped_data,
                indices=[999],
                backend="sequential",
                progress_bar=False,
            )

    def test_refit_loky_backend(self, initial_result, grouped_data):
        """Refit works with loky (multiprocessing) backend."""
        result2 = initial_result.fit(
            grouped_data,
            indices=[0, 1],
            backend="loky",
            n_jobs=2,
            wlmin=1.5,
            wlmax=5.0,
            method="leastsq",
            progress_bar=False,
        )
        assert len(result2) == len(initial_result)
        # Non-refitted groups copied
        assert result2["2"] is initial_result["2"]
        assert result2["3"] is initial_result["3"]
        # Refitted groups are new
        assert result2["0"] is not initial_result["0"]
        assert result2["1"] is not initial_result["1"]
        # Refitted groups have valid params
        assert result2["0"].params["thickness"].value > 0
        assert result2["1"].params["thickness"].value > 0

    def test_get_summary_dataframe(self, initial_result):
        """_get_summary_dataframe returns a proper DataFrame."""
        df = initial_result._get_summary_dataframe()
        assert isinstance(df, pd.DataFrame)
        assert "index" in df.columns
        assert "redchi" in df.columns
        assert "thickness" in df.columns
        assert "thickness_err" in df.columns
        assert len(df) == len(initial_result)

    def test_query_with_param_errors(self, initial_result, grouped_data):
        """Query can reference parameter error columns."""
        # This should not raise
        df = initial_result._get_summary_dataframe()
        # Verify thickness_err column exists for querying
        assert "thickness_err" in df.columns

        result2 = initial_result.fit(
            grouped_data,
            query="thickness_err < 10",
            backend="sequential",
            wlmin=1.5,
            wlmax=5.0,
            method="leastsq",
            progress_bar=False,
        )
        assert len(result2) == len(initial_result)

    def test_group_shape_preserved(self, initial_result, grouped_data):
        """Refitted result preserves group_shape from original."""
        result2 = initial_result.fit(
            grouped_data,
            indices=[0],
            backend="sequential",
            wlmin=1.5,
            wlmax=5.0,
            method="leastsq",
            progress_bar=False,
        )
        assert result2.group_shape == initial_result.group_shape

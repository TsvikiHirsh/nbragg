"""
Tests for Data.from_grouped_dataframes() and Data.from_grouped_transmission().
"""
import pytest
import pandas as pd
import numpy as np
from nbragg import Data


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_counts_df(n=40, seed=None):
    rng = np.random.default_rng(seed)
    tof = np.arange(100, 100 + n * 10, 10, dtype=float)
    counts = rng.poisson(1000, n).astype(float)
    return pd.DataFrame({"tof": tof, "counts": counts, "err": np.sqrt(counts)})


def _make_transmission_df(n=40, seed=None):
    rng = np.random.default_rng(seed)
    wl = np.linspace(0.5, 5.0, n)
    trans = rng.uniform(0.3, 0.9, n)
    err = trans * 0.02
    return pd.DataFrame({"wavelength": wl, "trans": trans, "err": err})


# ---------------------------------------------------------------------------
# from_grouped_dataframes – int indices (1-D)
# ---------------------------------------------------------------------------

class TestFromGroupedDataframesInt:
    def _signal_ob(self):
        signal = {i: _make_counts_df(seed=i) for i in range(3)}
        openbeam = {i: _make_counts_df(seed=i + 100) for i in range(3)}
        return signal, openbeam

    def test_is_grouped(self):
        signal, openbeam = self._signal_ob()
        data = Data.from_grouped_dataframes(signal, openbeam)
        assert data.is_grouped

    def test_indices(self):
        signal, openbeam = self._signal_ob()
        data = Data.from_grouped_dataframes(signal, openbeam)
        assert set(data.indices) == {"0", "1", "2"}

    def test_groups_have_transmission_columns(self):
        signal, openbeam = self._signal_ob()
        data = Data.from_grouped_dataframes(signal, openbeam)
        for idx in data.indices:
            df = data.groups[idx]
            assert set(df.columns) >= {"wavelength", "trans", "err"}

    def test_groups_signal_stored(self):
        signal, openbeam = self._signal_ob()
        data = Data.from_grouped_dataframes(signal, openbeam)
        for idx in data.indices:
            assert idx in data.groups_signal
            assert idx in data.groups_openbeam

    def test_default_table_set(self):
        signal, openbeam = self._signal_ob()
        data = Data.from_grouped_dataframes(signal, openbeam)
        assert data.table is not None

    def test_group_shape_1d(self):
        signal, openbeam = self._signal_ob()
        data = Data.from_grouped_dataframes(signal, openbeam)
        assert data.group_shape == (3,)


# ---------------------------------------------------------------------------
# from_grouped_dataframes – shared openbeam
# ---------------------------------------------------------------------------

class TestFromGroupedDataframesSharedOB:
    def test_shared_openbeam_dataframe(self):
        signal = {i: _make_counts_df(seed=i) for i in range(4)}
        shared_ob = _make_counts_df(seed=999)
        data = Data.from_grouped_dataframes(signal, shared_ob)
        assert data.is_grouped
        assert len(data.indices) == 4

    def test_transmission_computed_for_all_groups(self):
        signal = {i: _make_counts_df(seed=i) for i in range(4)}
        shared_ob = _make_counts_df(seed=999)
        data = Data.from_grouped_dataframes(signal, shared_ob)
        for idx in data.indices:
            df = data.groups[idx]
            assert len(df) > 0
            assert df["trans"].notna().all()


# ---------------------------------------------------------------------------
# from_grouped_dataframes – 2-D tuple indices
# ---------------------------------------------------------------------------

class TestFromGroupedDataframes2D:
    def _make_2d(self):
        coords = [(0, 0), (0, 1), (1, 0), (1, 1)]
        signal = {c: _make_counts_df(seed=i) for i, c in enumerate(coords)}
        openbeam = {c: _make_counts_df(seed=i + 50) for i, c in enumerate(coords)}
        return signal, openbeam, coords

    def test_is_grouped_2d(self):
        signal, openbeam, _ = self._make_2d()
        data = Data.from_grouped_dataframes(signal, openbeam)
        assert data.is_grouped

    def test_indices_2d_format(self):
        signal, openbeam, coords = self._make_2d()
        data = Data.from_grouped_dataframes(signal, openbeam)
        expected = {str(c).replace(" ", "") for c in coords}
        assert set(data.indices) == expected

    def test_group_shape_2d(self):
        signal, openbeam, _ = self._make_2d()
        data = Data.from_grouped_dataframes(signal, openbeam)
        assert data.group_shape is not None
        assert len(data.group_shape) == 2


# ---------------------------------------------------------------------------
# from_grouped_dataframes – string (named) indices
# ---------------------------------------------------------------------------

class TestFromGroupedDataframesNamed:
    def test_named_indices(self):
        names = ["sample_a", "sample_b", "reference"]
        signal = {n: _make_counts_df(seed=i) for i, n in enumerate(names)}
        openbeam = {n: _make_counts_df(seed=i + 20) for i, n in enumerate(names)}
        data = Data.from_grouped_dataframes(signal, openbeam)
        assert set(data.indices) == set(names)

    def test_group_shape_named_is_none(self):
        names = ["a", "b"]
        signal = {n: _make_counts_df(seed=i) for i, n in enumerate(names)}
        openbeam = {n: _make_counts_df(seed=i + 10) for i, n in enumerate(names)}
        data = Data.from_grouped_dataframes(signal, openbeam)
        assert data.group_shape is None


# ---------------------------------------------------------------------------
# from_grouped_dataframes – error handling
# ---------------------------------------------------------------------------

class TestFromGroupedDataframesErrors:
    def test_signal_not_dict_raises(self):
        with pytest.raises(TypeError, match="dict"):
            Data.from_grouped_dataframes([], {})

    def test_openbeam_wrong_type_raises(self):
        signal = {0: _make_counts_df()}
        with pytest.raises(TypeError, match="dict"):
            Data.from_grouped_dataframes(signal, "not_a_df_or_dict")

    def test_empty_signal_raises(self):
        with pytest.raises(ValueError, match="empty"):
            Data.from_grouped_dataframes({}, {})


# ---------------------------------------------------------------------------
# from_grouped_transmission
# ---------------------------------------------------------------------------

class TestFromGroupedTransmission:
    def _groups(self, n=3):
        return {i: _make_transmission_df(seed=i) for i in range(n)}

    def test_is_grouped(self):
        data = Data.from_grouped_transmission(self._groups())
        assert data.is_grouped

    def test_indices(self):
        data = Data.from_grouped_transmission(self._groups(4))
        assert set(data.indices) == {"0", "1", "2", "3"}

    def test_groups_have_correct_columns(self):
        data = Data.from_grouped_transmission(self._groups())
        for idx in data.indices:
            df = data.groups[idx]
            assert "wavelength" in df.columns
            assert "trans" in df.columns
            assert "err" in df.columns

    def test_default_table_set(self):
        data = Data.from_grouped_transmission(self._groups())
        assert data.table is not None

    def test_transmission_values_preserved(self):
        groups = self._groups(2)
        data = Data.from_grouped_transmission(groups)
        pd.testing.assert_frame_equal(
            data.groups["0"].reset_index(drop=True),
            groups[0].rename(columns={groups[0].columns[0]: "wavelength"}).reset_index(drop=True),
            check_like=True,
        )

    def test_2d_tuple_indices(self):
        coords = [(0, 0), (0, 1), (1, 0)]
        groups = {c: _make_transmission_df(seed=i) for i, c in enumerate(coords)}
        data = Data.from_grouped_transmission(groups)
        expected = {str(c).replace(" ", "") for c in coords}
        assert set(data.indices) == expected

    def test_named_indices(self):
        groups = {"roi_a": _make_transmission_df(seed=0), "roi_b": _make_transmission_df(seed=1)}
        data = Data.from_grouped_transmission(groups)
        assert set(data.indices) == {"roi_a", "roi_b"}

    def test_not_dict_raises(self):
        with pytest.raises(TypeError, match="dict"):
            Data.from_grouped_transmission([])

    def test_empty_dict_raises(self):
        with pytest.raises(ValueError, match="empty"):
            Data.from_grouped_transmission({})

    def test_dropna(self):
        groups = {0: _make_transmission_df(seed=0)}
        groups[0].loc[5, "trans"] = np.nan
        data = Data.from_grouped_transmission(groups, dropna=True)
        assert data.groups["0"]["trans"].notna().all()

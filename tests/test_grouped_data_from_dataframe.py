"""
Tests for Data.from_grouped with dict/list-of-DataFrame inputs.
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


# ---------------------------------------------------------------------------
# dict of DataFrames – int keys (1-D)
# ---------------------------------------------------------------------------

class TestDictInputInt:
    def _sig_ob(self):
        signal = {i: _make_counts_df(seed=i) for i in range(3)}
        openbeam = {i: _make_counts_df(seed=i + 100) for i in range(3)}
        return signal, openbeam

    def test_is_grouped(self):
        s, ob = self._sig_ob()
        data = Data.from_grouped(s, ob)
        assert data.is_grouped

    def test_indices(self):
        s, ob = self._sig_ob()
        data = Data.from_grouped(s, ob)
        assert set(data.indices) == {"0", "1", "2"}

    def test_groups_have_transmission_columns(self):
        s, ob = self._sig_ob()
        data = Data.from_grouped(s, ob)
        for idx in data.indices:
            df = data.groups[idx]
            assert {"wavelength", "trans", "err"}.issubset(df.columns)

    def test_groups_signal_stored(self):
        s, ob = self._sig_ob()
        data = Data.from_grouped(s, ob)
        for idx in data.indices:
            assert idx in data.groups_signal
            assert idx in data.groups_openbeam

    def test_default_table_set(self):
        s, ob = self._sig_ob()
        data = Data.from_grouped(s, ob)
        assert data.table is not None

    def test_group_shape_1d(self):
        s, ob = self._sig_ob()
        data = Data.from_grouped(s, ob)
        assert data.group_shape == (3,)


# ---------------------------------------------------------------------------
# dict of DataFrames – shared openbeam DataFrame
# ---------------------------------------------------------------------------

class TestDictInputSharedOB:
    def test_shared_openbeam(self):
        signal = {i: _make_counts_df(seed=i) for i in range(4)}
        shared_ob = _make_counts_df(seed=999)
        data = Data.from_grouped(signal, shared_ob)
        assert data.is_grouped
        assert len(data.indices) == 4

    def test_transmission_computed_for_all_groups(self):
        signal = {i: _make_counts_df(seed=i) for i in range(4)}
        shared_ob = _make_counts_df(seed=999)
        data = Data.from_grouped(signal, shared_ob)
        for idx in data.indices:
            assert data.groups[idx]["trans"].notna().all()


# ---------------------------------------------------------------------------
# dict of DataFrames – 2-D tuple keys
# ---------------------------------------------------------------------------

class TestDictInput2D:
    def _make_2d(self):
        coords = [(0, 0), (0, 1), (1, 0), (1, 1)]
        signal = {c: _make_counts_df(seed=i) for i, c in enumerate(coords)}
        openbeam = {c: _make_counts_df(seed=i + 50) for i, c in enumerate(coords)}
        return signal, openbeam, coords

    def test_indices_2d_format(self):
        s, ob, coords = self._make_2d()
        data = Data.from_grouped(s, ob)
        expected = {str(c).replace(" ", "") for c in coords}
        assert set(data.indices) == expected

    def test_group_shape_2d(self):
        s, ob, _ = self._make_2d()
        data = Data.from_grouped(s, ob)
        assert data.group_shape is not None
        assert len(data.group_shape) == 2


# ---------------------------------------------------------------------------
# dict of DataFrames – string (named) keys
# ---------------------------------------------------------------------------

class TestDictInputNamed:
    def test_named_indices(self):
        names = ["sample_a", "sample_b", "reference"]
        signal = {n: _make_counts_df(seed=i) for i, n in enumerate(names)}
        openbeam = {n: _make_counts_df(seed=i + 20) for i, n in enumerate(names)}
        data = Data.from_grouped(signal, openbeam)
        assert set(data.indices) == set(names)

    def test_group_shape_named_is_none(self):
        signal = {"a": _make_counts_df(seed=0), "b": _make_counts_df(seed=1)}
        openbeam = {"a": _make_counts_df(seed=10), "b": _make_counts_df(seed=11)}
        data = Data.from_grouped(signal, openbeam)
        assert data.group_shape is None


# ---------------------------------------------------------------------------
# dict of DataFrames – indices parameter overrides dict keys
# ---------------------------------------------------------------------------

class TestDictInputCustomIndices:
    def test_custom_indices_override_keys(self):
        signal = {0: _make_counts_df(seed=0), 1: _make_counts_df(seed=1)}
        openbeam = {0: _make_counts_df(seed=10), 1: _make_counts_df(seed=11)}
        data = Data.from_grouped(signal, openbeam, indices=["a", "b"])
        assert set(data.indices) == {"a", "b"}


# ---------------------------------------------------------------------------
# list of DataFrames
# ---------------------------------------------------------------------------

class TestListInput:
    def test_list_with_indices(self):
        sig_list = [_make_counts_df(seed=i) for i in range(3)]
        ob_list  = [_make_counts_df(seed=i + 10) for i in range(3)]
        data = Data.from_grouped(sig_list, ob_list, indices=["x", "y", "z"])
        assert data.is_grouped
        assert set(data.indices) == {"x", "y", "z"}

    def test_list_auto_int_indices(self):
        sig_list = [_make_counts_df(seed=i) for i in range(3)]
        ob_list  = [_make_counts_df(seed=i + 10) for i in range(3)]
        data = Data.from_grouped(sig_list, ob_list)
        assert set(data.indices) == {"0", "1", "2"}

    def test_list_shared_openbeam(self):
        sig_list = [_make_counts_df(seed=i) for i in range(3)]
        shared_ob = _make_counts_df(seed=99)
        data = Data.from_grouped(sig_list, shared_ob)
        assert len(data.indices) == 3


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

class TestErrorHandling:
    def test_empty_dict_raises(self):
        with pytest.raises(ValueError, match="empty"):
            Data.from_grouped({}, {})

    def test_openbeam_length_mismatch_raises(self):
        signal = {i: _make_counts_df(seed=i) for i in range(3)}
        openbeam = {i: _make_counts_df(seed=i) for i in range(2)}  # missing key 2
        with pytest.raises(ValueError):
            Data.from_grouped(signal, openbeam)

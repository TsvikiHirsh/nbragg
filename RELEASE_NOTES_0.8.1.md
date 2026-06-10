# nbragg v0.8.1 — Wl-resolution-invariant response + in-memory grouped-data constructors

This release adds two features:
1. A new response kind `jorgensen_inv` that produces a physical resolution width independent of how the user samples the wavelength grid.
2. In-memory pandas DataFrame inputs for `Data.from_grouped` and `Data.from_transmission`.

## Highlights

### `jorgensen_inv` — wavelength-resolution-invariant Jorgensen response

The existing `jorgensen` and `full_jorgensen` responses are convolved with the cross-section via `scipy.ndimage.convolve1d`, which is unit-agnostic: it applies the kernel sample-by-sample. As a result, the **physical** smearing width was `N_response_samples × data_dwl` — it scaled linearly with how finely the user sampled `wl`. Coarsen the data by 10× and the modelled "instrument resolution" widened by 10×.

The new `jorgensen_inv` kind evaluates the response on a Δλ grid whose step matches the actual data wavelength step (`data_dwl`). The physical kernel width is set by `α`, `β`, `σ` alone — independent of sampling.

```python
m = nbragg.TransmissionModel(
    xs,
    response="jorgensen_inv",
    vary_response=True,
)
m.params["α0"].set(value=23)
m.params["β0"].set(value=6)
```

Demonstration:

| samples in [3.6, 3.63] | `jorgensen` (old) | `jorgensen_inv` (new) |
|---|---|---|
| max\|T_coarse − T_fine_sub\| | **0.0194** | **0.00005** (~390× better) |

Test pinned in `tests/test_response.py::TestTransmissionModelInvariance`.

**Backwards-compatible.** Existing `jorgensen`/`full_jorgensen`/`square`/etc. behavior is unchanged. The `data_wl`/`data_dwl` kwargs are absorbed by the legacy responses via `**kwargs`. If you were tuning `α`/`β` against the old grid-dependent behavior, the converged values won't carry over to `jorgensen_inv`; expect to re-tune to your instrument's true resolution.

### `Data.from_grouped` now accepts dicts and lists of DataFrames
The classmethod that used to require glob patterns or folder paths now also accepts:

- **dict of DataFrames** — keys become the group indices
  ```python
  data = Data.from_grouped({0: sig0, 1: sig1}, {0: ob0, 1: ob1})
  ```

- **dict + shared openbeam** — pass a single DataFrame as `openbeam` and it is reused for every group
  ```python
  data = Data.from_grouped({0: sig0, 1: sig1}, shared_ob_df)
  ```

- **list of DataFrames** — sequential integer indices are assigned automatically; pass `indices=` to override
  ```python
  data = Data.from_grouped([sig0, sig1], [ob0, ob1], indices=["a", "b"])
  ```

- **2D / named indices** — the same `indices` semantics that the file-based path supports (`[(x, y), ...]` tuples for 2D grids, strings for named groups) work for in-memory input too.

`empty_signal` and `empty_openbeam` accept the same dict / list / single-DataFrame shapes.

### `Data.from_transmission` accepts dict / list of DataFrames
The transmission-side constructor mirrors the same input flexibility:
```python
data = Data.from_transmission({0: df0, 1: df1})
data = Data.from_transmission([df0, df1], indices=["roi_a", "roi_b"])
```

### Tests
- New `tests/test_grouped_data_from_dataframe.py` (28 tests) covers dict/list inputs, shared openbeam, 2D and named keys, custom `indices=`, error handling, and the same for `from_transmission`.
- New tests in `tests/test_response.py`: physical-width invariance across `data_dwl`, `data_wl` array handling, end-to-end T-invariance via `TransmissionModel`.

Full suite: **250/250 pass** (5 unrelated env-dependent skips).

## Backwards compatibility
Fully backwards-compatible. File-glob / folder-path inputs to `from_grouped` / `from_transmission` continue to work unchanged; legacy `jorgensen` and `full_jorgensen` responses are untouched.

## Upgrade
```bash
pip install --upgrade nbragg
```

### `Data.from_grouped` now accepts dicts and lists of DataFrames
The classmethod that used to require glob patterns or folder paths now also accepts:

- **dict of DataFrames** — keys become the group indices
  ```python
  data = Data.from_grouped({0: sig0, 1: sig1}, {0: ob0, 1: ob1})
  ```

- **dict + shared openbeam** — pass a single DataFrame as `openbeam` and it is reused for every group
  ```python
  data = Data.from_grouped({0: sig0, 1: sig1}, shared_ob_df)
  ```

- **list of DataFrames** — sequential integer indices are assigned automatically; pass `indices=` to override
  ```python
  data = Data.from_grouped([sig0, sig1], [ob0, ob1], indices=["a", "b"])
  ```

- **2D / named indices** — the same `indices` semantics that the file-based path supports (`[(x, y), ...]` tuples for 2D grids, strings for named groups) work for in-memory input too.

`empty_signal` and `empty_openbeam` accept the same dict / list / single-DataFrame shapes.

### `Data.from_transmission` accepts dict / list of DataFrames
The transmission-side constructor mirrors the same input flexibility:
```python
data = Data.from_transmission({0: df0, 1: df1})
data = Data.from_transmission([df0, df1], indices=["roi_a", "roi_b"])
```
Each DataFrame is expected to expose the same columns as the file-based path (`wavelength`/`tof`, `trans`, `err`).

### Tests
New `tests/test_grouped_data_from_dataframe.py` (28 tests) covers:
- Dict input with int keys, with a shared openbeam, with 2D-tuple keys, and with named string keys
- Custom `indices=` overriding dict keys
- List input with auto-assigned int indices and with explicit indices
- Shared openbeam DataFrame across a list of signal DataFrames
- Error handling: empty collections, openbeam-length mismatch
- `from_transmission` grouped variants matching the same shape rules

Full suite: **246/246 pass** (5 unrelated env-dependent skips).

## Backwards compatibility
Fully backwards-compatible. The file-glob and folder-path inputs to `from_grouped` / `from_transmission` continue to work unchanged — the new code path activates only when a dict, list, or single DataFrame is passed in.

## Upgrade
```bash
pip install --upgrade nbragg
```

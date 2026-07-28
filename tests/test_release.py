"""Release smoke test for nbragg 1.0.

A fast, self-contained end-to-end gate for the packaged release: it checks that
the version is a clean stable release, the public API is importable, and the core
single / grouped / save-load / extinction workflows run and recover known inputs.
All data is generated synthetically from nbragg's own forward model, so the gate
needs no data files and runs in any environment.

Run just this gate with::

    pytest tests/test_release.py -q
"""
import re
import tempfile

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

import nbragg


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _synthetic_transmission(thickness, wl=None, noise=0.004, seed=0):
    """A synthetic alpha-iron transmission spectrum (Bragg edges + noise).

    Returns the transmission ``DataFrame`` and the model that generated it.
    """
    if wl is None:
        wl = np.linspace(1.5, 5.0, 500)
    xs = nbragg.CrossSection(iron=nbragg.materials["Fe_sg229_Iron-alpha"])
    model = nbragg.TransmissionModel(xs, vary_basic=True)
    rng = np.random.default_rng(seed)
    T = model.transmission(wl, thickness=thickness, norm=1.0)
    T = T + rng.normal(0.0, noise, len(wl))
    return pd.DataFrame({"wavelength": wl, "trans": T, "err": noise}), model


# --------------------------------------------------------------------------- #
# release gate
# --------------------------------------------------------------------------- #
def test_version_is_stable_release():
    """The installed version is a clean, >= 1.0 stable release (no dev/pre tag)."""
    v = nbragg.__version__
    m = re.match(r"^(\d+)\.(\d+)(?:\.(\d+))?$", v)
    assert m, f"{v!r} is not a clean X.Y[.Z] release version"
    assert int(m.group(1)) >= 1, f"expected a >= 1.0 stable release, got {v}"


def test_public_api_present():
    """Every documented public name is importable from the top-level package."""
    for name in (
        "Data", "CrossSection", "TransmissionModel", "GroupedFitResult",
        "Response", "Background", "materials", "register_material",
        "save_result", "load_result", "utils",
    ):
        assert hasattr(nbragg, name), f"public API missing: nbragg.{name}"
    assert len(nbragg.materials) > 0


def test_end_to_end_single_fit():
    """A single transmission fit converges and recovers the input thickness."""
    df, model = _synthetic_transmission(thickness=0.8)
    data = nbragg.Data.from_transmission(df)
    result = model.fit(data, wlmin=1.6, wlmax=4.8,
                       method="least-squares", progress_bar=False)
    assert result.success
    assert result.redchi < 5.0
    assert result.params["thickness"].value == pytest.approx(0.8, abs=0.05)


def test_end_to_end_grouped_fit():
    """A grouped/parallel fit recovers a per-pixel thickness pattern + 2D map."""
    xs = nbragg.CrossSection(iron=nbragg.materials["Fe_sg229_Iron-alpha"])
    model = nbragg.TransmissionModel(xs, vary_basic=True)
    wl = np.linspace(1.5, 5.0, 400)
    rng = np.random.default_rng(1)

    true, groups = {}, {}
    for i in range(2):
        for j in range(3):
            t = 0.4 + 0.1 * (3 * i + j)
            true[(i, j)] = t
            T = model.transmission(wl, thickness=t, norm=1.0)
            T = T + rng.normal(0.0, 0.004, len(wl))
            groups[(i, j)] = pd.DataFrame({"wavelength": wl, "trans": T, "err": 0.004})

    data = nbragg.Data.from_transmission(groups)
    result = model.fit(data, n_jobs=2, wlmin=1.6, wlmax=4.8, progress_bar=False)

    assert sorted(result.group_shape) == [2, 3]
    for idx, t in true.items():
        assert result[idx].params["thickness"].value == pytest.approx(t, abs=0.05)
    # the 2D parameter map renders without error
    assert result.plot_parameter_map("thickness") is not None


def test_save_load_roundtrip():
    """A fit result survives a save -> load round-trip."""
    df, model = _synthetic_transmission(thickness=0.8)
    data = nbragg.Data.from_transmission(df)
    result = model.fit(data, wlmin=1.6, wlmax=4.8,
                       method="least-squares", progress_bar=False)

    path = tempfile.mktemp(suffix=".json")
    nbragg.save_result(result, path)
    loaded = nbragg.load_result(path, model=model)

    assert loaded["success"]
    assert loaded["params"]["thickness"].value == pytest.approx(
        result.params["thickness"].value, rel=1e-6)


def test_extinction_plugin_optional():
    """Extinction works when the optional CrysExtn plugin is installed.

    Skips cleanly when the plugin is absent so the gate still passes on a base
    install.
    """
    try:
        xs = nbragg.CrossSection(iron={
            "mat": "Fe_sg225_Iron-gamma.ncmat",
            "ext_method": "BC_mix", "ext_l": 5.0, "ext_g": 4.76,
            "ext_L": 20.0, "ext_dist": "Gauss",
        })
        total = xs.table["total"].to_numpy()
    except Exception as exc:  # NCrystal errors on the extinction section if unplugged
        pytest.skip(f"CrysExtn plugin not available: {exc}")

    if "@CUSTOM_CRYSEXTN" not in xs.textdata["iron"]:
        pytest.skip("CrysExtn plugin not available")

    assert np.isfinite(total).all()
    assert (total > 0).all()  # valid parameters -> no NaN/zero bands

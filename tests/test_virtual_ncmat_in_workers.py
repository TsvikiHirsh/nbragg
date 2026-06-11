"""
Tests for the in-memory NCrystal data fix for grouped fitting.

Issue: NCrystal's in-memory file registry (NC.registerInMemoryFileData) lives
in the C++ singleton per-process. When grouped fitting spawns worker processes,
those workers see an empty virtual-file registry — any model that references
a virtual ncmat (either user-registered or one nbragg created internally for
lattice/extinction updates) fails with NCFileNotFound.

Fix: snapshot the virtual registry in the parent before spawning workers, then
re-register everything in each worker via a ProcessPoolExecutor initializer.
"""
import multiprocessing as mp
import os

import numpy as np
import pandas as pd
import pytest

import NCrystal as NC


def _read_text(path):
    with open(path, "r") as f:
        return f.read()


@pytest.fixture
def iron_alpha_content():
    """Read a real ncmat file off disk so we can also register it as virtual."""
    path = os.path.join(os.path.dirname(__file__), "Fe_sg229_Iron-alpha_CrysExtn1.ncmat")
    return _read_text(path)


def test_snapshot_and_re_register_roundtrip(iron_alpha_content):
    """The two helpers must round-trip virtual content faithfully."""
    from nbragg.grouped_fit import _snapshot_virtual_ncmat_files, _init_worker_virtuals

    NC.registerInMemoryFileData("test_virt_iron.ncmat", iron_alpha_content)

    snapshot = _snapshot_virtual_ncmat_files()
    assert "test_virt_iron.ncmat" in snapshot
    assert snapshot["test_virt_iron.ncmat"].startswith(iron_alpha_content[:20])

    # Re-registering (e.g. in a worker) is a no-op for content equality
    _init_worker_virtuals(snapshot)
    td = NC.createTextData("virtual::test_virt_iron.ncmat")
    assert td.rawData == iron_alpha_content


def _worker_check_virtual(name):
    """Run in a child process: try to load the named virtual file."""
    import NCrystal as NC

    try:
        td = NC.createTextData(f"virtual::{name}")
        return ("ok", len(td.rawData))
    except Exception as e:
        return ("err", f"{type(e).__name__}: {e}")


def test_initializer_makes_virtual_files_visible_in_child(iron_alpha_content):
    """The actual end-to-end fix: with the initializer, a child process can
    resolve a virtual file the parent registered. Without it, it can't."""
    from concurrent.futures import ProcessPoolExecutor
    from nbragg.grouped_fit import _snapshot_virtual_ncmat_files, _init_worker_virtuals

    NC.registerInMemoryFileData("xfer_check.ncmat", iron_alpha_content)
    snapshot = _snapshot_virtual_ncmat_files()
    assert "xfer_check.ncmat" in snapshot

    # WITHOUT the initializer: the child sees no virtuals (NCFileNotFound).
    # Use 'spawn' to guarantee a fresh interpreter — most reliable on Linux.
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=1, mp_context=ctx) as ex:
        status, _ = ex.submit(_worker_check_virtual, "xfer_check.ncmat").result(timeout=30)
    assert status == "err", "control case: child should NOT see virtual without initializer"

    # WITH the initializer: child resolves it.
    with ProcessPoolExecutor(
        max_workers=1,
        mp_context=ctx,
        initializer=_init_worker_virtuals,
        initargs=(snapshot,),
    ) as ex:
        status, n = ex.submit(_worker_check_virtual, "xfer_check.ncmat").result(timeout=30)
    assert status == "ok", "with initializer: child must see virtual file"
    assert n == len(iron_alpha_content)


def test_empty_snapshot_is_harmless():
    """An empty snapshot must not crash the initializer or workers."""
    from nbragg.grouped_fit import _init_worker_virtuals
    _init_worker_virtuals({})       # explicit empty
    _init_worker_virtuals(None)     # tolerate falsy

# nbragg v0.8.2 — Grouped fits work with in-memory NCrystal data

Patch release that makes parallel grouped fitting (`_fit_grouped_loky`, `_refit_loky`) work correctly when the model references in-memory NCrystal materials — either user-registered via `NC.registerInMemoryFileData`, or virtual `.nbragg` template files that nbragg itself creates internally for lattice/extinction edits.

## The bug

Grouped fits spawn workers via `concurrent.futures.ProcessPoolExecutor`. NCrystal's in-memory file registry lives in a C++ singleton **per process** — it does not survive fork/spawn. Workers therefore booted with an empty virtual-file registry, and any model that referenced a virtual ncmat failed with `NCFileNotFound` (e.g. `'Could not find data: "my_virt.ncmat"'`).

The visible symptom: as soon as users either (a) registered their own ncmat content in memory before fitting, or (b) used certain extinction/lattice features that cause nbragg to mint a virtual template file, parallel grouped fits failed for every group.

## The fix

Before launching workers, the parent process now snapshots every virtual file in NCrystal's registry using `NC.browseFiles(factory="virtual")` plus `NC.createTextData(entry.fullKey).rawData`. The snapshot is passed to `ProcessPoolExecutor(initializer=..., initargs=(snapshot,))`; the initializer re-registers each file in the worker before it touches the model.

- [src/nbragg/grouped_fit.py](src/nbragg/grouped_fit.py) — new module-level helpers `_snapshot_virtual_ncmat_files` and `_init_worker_virtuals`, wired into `_refit_loky`.
- [src/nbragg/model_fitting.py](src/nbragg/model_fitting.py) — same wiring for `_fit_grouped_loky`.

No changes to user-facing APIs. Empty snapshots are a no-op, so the fix adds no cost when no virtual files are registered.

## Tests
Three new tests in [tests/test_virtual_ncmat_in_workers.py](tests/test_virtual_ncmat_in_workers.py):
- **Snapshot/re-register roundtrip** preserves content byte-for-byte.
- **End-to-end control:** with a fresh `mp.get_context("spawn")` pool, the child process **cannot** see the parent's virtual file *without* the initializer (asserted), and **can** see it *with* the initializer (asserted). This pins the fix to the exact mechanism that was broken.
- **Empty-snapshot no-op** path.

Full suite: **253/253 pass** (was 250 in 0.8.1).

## Backwards compatibility
Fully backwards-compatible. No API changes; behavior only differs for code paths that were previously broken (where the model referenced a virtual ncmat).

## Upgrade
```bash
pip install --upgrade nbragg
```

## Credits
Reported by [@wolfertz](https://github.com/wolfertz) with a clear reproduction and the correct intuition that NCrystal's process-local registry was the culprit.

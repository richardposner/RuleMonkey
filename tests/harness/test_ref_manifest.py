#!/usr/bin/env python3
"""Checks on `harness/_ref_manifest.py`, the reference-tree integrity gate.

Two things are under test here.  The unit half pins what a MANIFEST.tsv
covers: the vendored reference artifacts, and not the gitignored
`replicates/` scratch they were aggregated from — a manifest that hashes
the scratch cannot verify on a clean clone, which is what took the
basicmodels parity suite out of reach (issue #63).

The other half verifies the manifests actually committed to this
repository against the live tree, from whatever checkout is running.  On
CI that checkout is a clean clone, so this is the guard that keeps every
parity harness startable there.

Pure stdlib, no rm_driver: runs on every ctest leg.

Usage:
    python3 tests/harness/test_ref_manifest.py
"""

from __future__ import annotations

import contextlib
import io
import os
import shutil
import sys
import tempfile

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "harness"))

import _ref_manifest as rm  # noqa: E402


def write(path: str, text: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def make_tree(root: str, *, with_scratch: bool) -> None:
    """A miniature reference root: vendored aggregates, optional scratch."""
    write(os.path.join(root, "PROVENANCE.md"), "# provenance\n")
    write(os.path.join(root, "sim_params.tsv"), "model\tt_end\nr01\t10\n")
    write(os.path.join(root, "ensemble", "r01.mean.tsv"), "time\tA\n0\t100\n")
    write(os.path.join(root, "ensemble", "r01.std.tsv"), "time\tA\n0\t0\n")
    if with_scratch:
        for rep in range(1, 4):
            write(os.path.join(root, "replicates", "r01", f"rep_{rep:03d}.gdat"), f"# {rep}\n")


VENDORED = {
    "PROVENANCE.md",
    "sim_params.tsv",
    "ensemble/r01.mean.tsv",
    "ensemble/r01.std.tsv",
}


def test_manifest_lists_vendored_files_only() -> None:
    """`replicates/` and OS noise stay out of a freshly written manifest."""
    with tempfile.TemporaryDirectory() as tmp:
        make_tree(tmp, with_scratch=True)
        write(os.path.join(tmp, ".DS_Store"), "mac noise\n")
        write(os.path.join(tmp, "ensemble", ".DS_Store"), "mac noise\n")
        entries = rm.read_manifest(rm.write_manifest(tmp))
        assert set(entries) == VENDORED, f"manifest covers {sorted(entries)}"


def test_manifest_paths_are_slash_separated() -> None:
    """A manifest written on Windows has to verify on Linux and vice versa."""
    with tempfile.TemporaryDirectory() as tmp:
        make_tree(tmp, with_scratch=False)
        with open(rm.write_manifest(tmp)) as f:
            body = [ln for ln in f.read().splitlines() if ln and not ln.startswith("#")]
        nested = [ln for ln in body if ln.startswith("ensemble")]
        assert len(nested) == 2, f"expected 2 nested rows, got {body}"
        assert not any("\\" in ln.split("\t")[0] for ln in body), body


def test_verify_passes_with_and_without_scratch() -> None:
    """The same manifest verifies on a clean clone and after a regen run."""
    with tempfile.TemporaryDirectory() as tmp:
        make_tree(tmp, with_scratch=True)
        rm.write_manifest(tmp)
        ok, problems = rm.verify_manifest(tmp)
        assert ok, f"scratch present: {problems}"

        shutil.rmtree(os.path.join(tmp, "replicates"))
        ok, problems = rm.verify_manifest(tmp)
        assert ok, f"scratch absent: {problems}"

        # ... and appears again once someone regenerates the ensemble.
        make_tree(tmp, with_scratch=True)
        ok, problems = rm.verify_manifest(tmp)
        assert ok, f"scratch regenerated: {problems}"


def reported(problems: list[str], needle: str) -> bool:
    return any(needle in p for p in problems)


def test_verify_catches_drift_in_vendored_files() -> None:
    """The point of the gate: an edited, missing or added reference is caught."""
    with tempfile.TemporaryDirectory() as tmp:
        make_tree(tmp, with_scratch=False)
        rm.write_manifest(tmp)
        mean = os.path.join(tmp, "ensemble", "r01.mean.tsv")

        write(mean, "time\tA\n0\t101\n")
        ok, problems = rm.verify_manifest(tmp)
        assert not ok and reported(problems, "hash mismatch for ensemble/r01.mean"), problems

        os.remove(mean)
        ok, problems = rm.verify_manifest(tmp)
        assert not ok and reported(problems, "missing reference file: ensemble/r01"), problems

        write(os.path.join(tmp, "ensemble", "r02.mean.tsv"), "time\tA\n0\t7\n")
        ok, problems = rm.verify_manifest(tmp)
        assert not ok and reported(problems, "untracked file in reference tree"), problems


def test_manifest_carrying_scratch_rows_is_one_diagnostic() -> None:
    """The issue-#63 shape: stale rows say so, once, instead of 2900 times."""
    with tempfile.TemporaryDirectory() as tmp:
        make_tree(tmp, with_scratch=False)
        manifest = rm.write_manifest(tmp)
        with open(manifest, "a") as f:
            for rep in range(1, 101):
                f.write(f"replicates/r01/rep_{rep:03d}.gdat\t{'0' * 64}\n")

        ok, problems = rm.verify_manifest(tmp)
        assert not ok, "a manifest reaching outside its coverage should not verify"
        assert len(problems) == 1, problems
        assert "100 path(s) outside what it covers" in problems[0], problems[0]
        assert "replicates/r01/rep_001.gdat" in problems[0], problems[0]


def test_enforce_exits_2_with_the_callers_own_regen_flags() -> None:
    """The abort names flags the calling script actually accepts."""
    with tempfile.TemporaryDirectory() as tmp:
        make_tree(tmp, with_scratch=False)
        rm.write_manifest(tmp)
        os.remove(os.path.join(tmp, "sim_params.tsv"))

        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            try:
                rm.enforce_or_warn(
                    tmp,
                    strict=True,
                    label="unit_test",
                    regen_hint="re-run with --no-verify-manifest --write-manifest",
                )
            except SystemExit as exc:
                code = exc.code
            else:
                raise AssertionError("strict enforce_or_warn returned instead of exiting")
        assert code == 2, code
        text = err.getvalue()
        assert "--no-verify-manifest --write-manifest" in text, text
        assert "--generate-refs" not in text, text

        # Warn mode reports the same problems and lets the caller continue.
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rm.enforce_or_warn(tmp, strict=False, label="unit_test")
        assert "missing reference file: sim_params.tsv" in err.getvalue(), err.getvalue()


def test_problem_list_is_capped() -> None:
    """A tree that lost a directory prints a summary, not one line per file."""
    with tempfile.TemporaryDirectory() as tmp:
        make_tree(tmp, with_scratch=False)
        for n in range(2, 60):
            write(os.path.join(tmp, "ensemble", f"r{n:02d}.mean.tsv"), f"time\tA\n0\t{n}\n")
        rm.write_manifest(tmp)
        shutil.rmtree(os.path.join(tmp, "ensemble"))

        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.suppress(SystemExit):
            rm.enforce_or_warn(tmp, strict=True, label="unit_test")
        bullets = [ln for ln in err.getvalue().splitlines() if ln.startswith("  - ")]
        assert len(bullets) == rm.MAX_REPORTED_PROBLEMS + 1, err.getvalue()
        assert bullets[-1].endswith("more"), bullets[-1]


def test_committed_reference_trees_verify() -> None:
    """Every MANIFEST.tsv in this checkout matches the tree next to it.

    On a clean clone this is the issue-#63 regression test: a manifest
    that gates on gitignored files fails here before anyone waits on a
    parity run to tell them.
    """
    roots = []
    for dirpath, dirnames, filenames in os.walk(os.path.join(REPO_ROOT, "tests")):
        dirnames[:] = [d for d in dirnames if d not in rm.EXCLUDED_DIRS]
        if rm.MANIFEST_FILENAME in filenames:
            roots.append(dirpath)
    assert roots, "no MANIFEST.tsv found under tests/"

    for root in sorted(roots):
        ok, problems = rm.verify_manifest(root)
        rel = os.path.relpath(root, REPO_ROOT)
        assert ok, f"{rel}: {len(problems)} problem(s); first: {problems[0]}"


TESTS = [
    test_manifest_lists_vendored_files_only,
    test_manifest_paths_are_slash_separated,
    test_verify_passes_with_and_without_scratch,
    test_verify_catches_drift_in_vendored_files,
    test_manifest_carrying_scratch_rows_is_one_diagnostic,
    test_enforce_exits_2_with_the_callers_own_regen_flags,
    test_problem_list_is_capped,
    test_committed_reference_trees_verify,
]


def main() -> int:
    failures = 0
    for test in TESTS:
        try:
            test()
        except AssertionError as exc:
            failures += 1
            print(f"FAIL  {test.__name__}\n      {exc}", flush=True)
        else:
            print(f"PASS  {test.__name__}", flush=True)
    print(f"TOTAL: {len(TESTS) - failures} PASS / {failures} FAIL", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

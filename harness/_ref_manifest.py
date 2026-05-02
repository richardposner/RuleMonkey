"""Reference-data manifest helpers shared by feature-coverage and basicmodels harnesses.

A MANIFEST.tsv records SHA-256 hashes of every reference file under a
given root directory.  The validate path verifies the live tree against
the manifest before running RM; mismatches abort with a clear diagnostic
so a stale or hand-edited reference can't silently change the verdict.

Format (tab-separated, sorted by relative path):

    # rulemonkey reference manifest
    # generated <UTC ISO 8601>
    # root <relative path of reference root from REPO_ROOT>
    <relative path>\t<sha256 hex>

Comment lines (`#`) are ignored on read.

Generation is opt-in (only fired by --generate-refs / --force-refs);
verification is the default at script start so a working tree where
references drifted gets caught.
"""

from __future__ import annotations

import datetime
import hashlib
import os
import sys

MANIFEST_FILENAME = "MANIFEST.tsv"


def _walk_files(root: str) -> list[str]:
    out = []
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            if fn == MANIFEST_FILENAME:
                continue
            out.append(os.path.join(dirpath, fn))
    out.sort()
    return out


def _hash_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(ref_root: str) -> str:
    """Write MANIFEST.tsv at `ref_root`.  Returns the manifest path."""
    if not os.path.isdir(ref_root):
        raise RuntimeError(f"reference root does not exist: {ref_root}")
    manifest_path = os.path.join(ref_root, MANIFEST_FILENAME)
    files = _walk_files(ref_root)
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(manifest_path, "w") as f:
        f.write("# rulemonkey reference manifest\n")
        f.write(f"# generated {now}\n")
        f.write(f"# root {os.path.relpath(ref_root)}\n")
        for path in files:
            rel = os.path.relpath(path, ref_root)
            f.write(f"{rel}\t{_hash_file(path)}\n")
    return manifest_path


def read_manifest(manifest_path: str) -> dict[str, str]:
    """Parse MANIFEST.tsv → dict of relative path → sha256 hex."""
    out = {}
    with open(manifest_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            out[parts[0]] = parts[1]
    return out


def verify_manifest(ref_root: str, *, strict: bool = True) -> tuple[bool, list[str]]:
    """Compare the live tree under `ref_root` against MANIFEST.tsv.

    Returns (ok, problems) where `problems` is a list of human-readable
    diagnostic strings.  When `strict` is False (warn-only mode) the
    caller may proceed despite mismatches; when True (default) the caller
    should refuse to run.
    """
    manifest_path = os.path.join(ref_root, MANIFEST_FILENAME)
    if not os.path.exists(manifest_path):
        return False, [
            f"no MANIFEST.tsv at {ref_root}; reference integrity is unverified — "
            f"regenerate refs with --generate-refs/--force-refs to bootstrap the manifest"
        ]
    expected = read_manifest(manifest_path)
    live_paths = {os.path.relpath(p, ref_root): p for p in _walk_files(ref_root)}

    problems = []
    for rel, want_hash in expected.items():
        abs_path = os.path.join(ref_root, rel)
        if not os.path.exists(abs_path):
            problems.append(f"missing reference file: {rel}")
            continue
        got = _hash_file(abs_path)
        if got != want_hash:
            problems.append(f"hash mismatch for {rel}: manifest={want_hash[:12]}… live={got[:12]}…")
    extras = sorted(set(live_paths) - set(expected))
    for rel in extras:
        problems.append(f"untracked file in reference tree: {rel} (not in MANIFEST.tsv)")

    return (not problems), problems


def enforce_or_warn(ref_root: str, *, strict: bool, label: str) -> None:
    """Verify the manifest at `ref_root` and act on the result.

    `strict=True` aborts the process on mismatch; `False` warns to stderr
    and returns.  Suitable as a startup gate for harness scripts.
    """
    ok, problems = verify_manifest(ref_root, strict=strict)
    if ok:
        return
    header = f"{label}: reference manifest verification reported {len(problems)} issue(s):"
    print(header, file=sys.stderr)
    for p in problems:
        print(f"  - {p}", file=sys.stderr)
    if strict:
        print(
            f"{label}: refusing to run validate path against an unverified reference tree.\n"
            f"To regenerate the manifest after intentional changes, re-run with "
            f"--generate-refs (or --force-refs) so the new state is committed.",
            file=sys.stderr,
        )
        sys.exit(2)

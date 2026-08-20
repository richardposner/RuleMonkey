"""Reference-data manifest helpers shared by feature-coverage and basicmodels harnesses.

A MANIFEST.tsv records SHA-256 hashes of the vendored reference files
under a given root directory (see "What the manifest covers" below for
what that excludes).  The validate path verifies the live tree against
the manifest before running RM; mismatches abort with a clear diagnostic
so a stale or hand-edited reference can't silently change the verdict.

Format (tab-separated, sorted by relative path):

    # rulemonkey reference manifest
    # generated <UTC ISO 8601>
    # root <path of the reference root, relative to the repo root>
    <relative path>\t<sha256 hex>

Relative paths are `/`-separated on every platform, so a manifest
written on one OS verifies on another.  Comment lines (`#`) are ignored
on read.

## What the manifest covers

The vendored reference artifacts, and only those.  Paths `.gitignore`
keeps out of the repository — the per-replicate ensembles under
`replicates/`, OS noise like `.DS_Store` — are skipped on both write and
verify.  A clean clone has none of those files and a machine that just
regenerated an ensemble has thousands of them, so hashing them makes the
gate fail one way in the first case ("missing reference file:
replicates/r01/rep_001.gdat", issue #63) and the other way in the second
("untracked file in reference tree").  Neither is reference drift: the
verdict path reads the aggregated `ensemble/*.{mean,std,tint}.tsv`, and
`replicates/` is the disposable scratch those were aggregated from.

Generation is opt-in (only fired by --generate-refs / --force-refs /
--write-manifest); verification is the default at script start so a
working tree where references drifted gets caught.
"""

from __future__ import annotations

import datetime
import hashlib
import os
import sys

MANIFEST_FILENAME = "MANIFEST.tsv"

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Mirrors `.gitignore` (`tests/reference/*/replicates/`, plus the OS-noise
# entries) — see "What the manifest covers" above for why these are skipped
# rather than hashed.
EXCLUDED_DIRS = frozenset({"replicates"})
EXCLUDED_FILES = frozenset({MANIFEST_FILENAME, ".DS_Store", "Thumbs.db"})

# Cap on the diagnostic list.  A ref tree missing one whole directory has
# one problem per file, and a few thousand of those scroll the closing
# line — the one saying what to do about it — off the screen.
MAX_REPORTED_PROBLEMS = 20


def _rel(path: str, root: str) -> str:
    """Manifest-relative path for `path`, `/`-separated on every platform."""
    return os.path.relpath(path, root).replace(os.sep, "/")


def _header_root(ref_root: str) -> str:
    """Value for the `# root` header line: repo-relative where it can be.

    Cosmetic — `read_manifest` skips comment lines — but it used to be
    `os.path.relpath(ref_root)`, relative to the *working directory*,
    which made the line depend on where the harness was invoked from and
    made it a hard error on Windows for any tree on another drive than
    the CWD (`ValueError: path is on mount 'C:', start on mount 'D:'`).
    A reference root outside the repository gets its absolute path.
    """
    ref_root = os.path.abspath(ref_root)
    try:
        rel = os.path.relpath(ref_root, REPO_ROOT)
    except ValueError:  # different Windows drive; no relative path exists
        rel = os.pardir
    if rel.split(os.sep)[0] == os.pardir:
        return ref_root.replace(os.sep, "/")
    return rel.replace(os.sep, "/")


def _is_excluded(rel: str) -> bool:
    """True for a manifest-relative path the manifest deliberately skips."""
    parts = rel.split("/")
    return parts[-1] in EXCLUDED_FILES or any(p in EXCLUDED_DIRS for p in parts[:-1])


def _walk_files(root: str) -> list[str]:
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in EXCLUDED_DIRS]
        for fn in filenames:
            if fn in EXCLUDED_FILES:
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
        f.write(f"# root {_header_root(ref_root)}\n")
        for path in files:
            f.write(f"{_rel(path, ref_root)}\t{_hash_file(path)}\n")
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
            f"no MANIFEST.tsv at {ref_root}; reference integrity is unverified until "
            f"the manifest is bootstrapped"
        ]
    expected = read_manifest(manifest_path)
    live_paths = {_rel(p, ref_root): p for p in _walk_files(ref_root)}

    problems = []

    # Rows the manifest should never have carried in the first place.
    # Reported once, not once per file: a manifest bootstrapped on a
    # machine that still had its `replicates/` scratch picks up thousands
    # of these, and each one is the same single fact about the manifest.
    uncovered = sorted(rel for rel in expected if _is_excluded(rel))
    if uncovered:
        problems.append(
            f"manifest lists {len(uncovered)} path(s) outside what it covers "
            f"(first: {uncovered[0]}); those are regenerated locally and gitignored — "
            f"rewrite the manifest so it lists the vendored files only"
        )

    for rel, want_hash in expected.items():
        if _is_excluded(rel):
            continue
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


def enforce_or_warn(ref_root: str, *, strict: bool, label: str, regen_hint: str = "") -> None:
    """Verify the manifest at `ref_root` and act on the result.

    `strict=True` aborts the process on mismatch; `False` warns to stderr
    and returns.  Suitable as a startup gate for harness scripts.
    `regen_hint` names the flags *this caller* accepts for rewriting the
    manifest, since they differ per script.
    """
    ok, problems = verify_manifest(ref_root, strict=strict)
    if ok:
        return
    header = f"{label}: reference manifest verification reported {len(problems)} issue(s):"
    print(header, file=sys.stderr)
    for p in problems[:MAX_REPORTED_PROBLEMS]:
        print(f"  - {p}", file=sys.stderr)
    if len(problems) > MAX_REPORTED_PROBLEMS:
        print(f"  - … and {len(problems) - MAX_REPORTED_PROBLEMS} more", file=sys.stderr)
    if strict:
        hint = regen_hint or "re-run the harness's reference-regeneration path"
        print(
            f"{label}: refusing to run validate path against an unverified reference tree.\n"
            f"To regenerate the manifest after intentional changes, {hint} "
            f"so the new state is committed.",
            file=sys.stderr,
        )
        sys.exit(2)

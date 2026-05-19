#!/usr/bin/env python3
"""Refresh, summarize, or drift-check the vendored BNGsim expression layer.

RuleMonkey's standalone rate-law / function / parameter expression evaluator
is BNGsim's ExprTkEvaluator (issue #6).  For non-BNGsim builds, RuleMonkey
vendors BNGsim's expression wrapper plus the stock ExprTk header that BNGsim
pins.  Inside a BNGsim build, CMake links the host bngsim::expression target
instead.

Refresh writes the vendored files and machine-readable provenance:

  scripts/vendor_exprtk.py --bngsim-repo /path/to/bngsim
  scripts/vendor_exprtk.py --bngsim-repo /path/to/bngsim --ref <commit>

No-write modes:

  scripts/vendor_exprtk.py --bngsim-repo /path/to/bngsim --summary
  scripts/vendor_exprtk.py --check
  scripts/vendor_exprtk.py --bngsim-repo /path/to/bngsim --check

Without --bngsim-repo, --check verifies the checked-in files against
third_party/bngsim_expr/VENDOR.json.  With --bngsim-repo, it also verifies
those files against the pinned BNGsim commit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


RM_ROOT = Path(__file__).resolve().parents[1]
VENDOR_DIR = RM_ROOT / "third_party" / "bngsim_expr"
VENDOR_NOTE = VENDOR_DIR / "VENDOR"
VENDOR_JSON = VENDOR_DIR / "VENDOR.json"
PIN_PREFIX = "Pinned commit: "
DEFAULT_REF = "HEAD"

# vendored path (RuleMonkey-relative) -> source path (BNGsim-checkout-relative)
FILES = {
    "exprtk_header": {
        "path": "third_party/exprtk/exprtk.hpp",
        "source_path": "third_party/exprtk/exprtk.hpp",
        "description": "stock ExprTk header pinned by BNGsim",
    },
    "expression_header": {
        "path": "third_party/bngsim_expr/include/bngsim/expression.hpp",
        "source_path": "include/bngsim/expression.hpp",
        "description": "BNGsim ExprTkEvaluator public wrapper header",
    },
    "expression_source": {
        "path": "third_party/bngsim_expr/src/expression.cpp",
        "source_path": "src/expression.cpp",
        "description": "BNGsim ExprTkEvaluator implementation",
    },
}


def run(cmd: list[str], *, text: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )


def git(bngsim: Path, args: list[str]) -> str:
    return run(["git", "-C", str(bngsim), *args]).stdout.strip()


def git_lines(bngsim: Path, args: list[str]) -> list[str]:
    output = git(bngsim, args)
    return [line for line in output.splitlines() if line]


def git_show(bngsim: Path, ref: str, rel: str) -> bytes:
    """Return bytes of a BNGsim-relative path at ref."""
    return run(["git", "-C", str(bngsim), "show", f"{ref}:./{rel}"], text=False).stdout


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def json_loads_bytes(data: bytes) -> dict[str, Any] | None:
    try:
        return json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return None


def remote_url(bngsim: Path, remote: str) -> str | None:
    try:
        return git(bngsim, ["remote", "get-url", remote])
    except subprocess.CalledProcessError:
        return None


def read_vendor_metadata() -> dict[str, Any] | None:
    if not VENDOR_JSON.exists():
        return None
    return json.loads(VENDOR_JSON.read_text())


def read_legacy_pin() -> str | None:
    if not VENDOR_NOTE.exists():
        return None
    for line in VENDOR_NOTE.read_text().splitlines():
        if line.startswith(PIN_PREFIX):
            return line[len(PIN_PREFIX) :].strip()
    return None


def resolve_bngsim_repo(path: Path | None, *, required: bool) -> Path | None:
    raw = path or (Path(os.environ["BNGSIM_REPO"]) if "BNGSIM_REPO" in os.environ else None)
    if raw is None:
        if required:
            raise RuntimeError("pass --bngsim-repo /path/to/bngsim or set BNGSIM_REPO")
        return None

    bngsim = raw.expanduser().resolve()
    if not bngsim.is_dir():
        raise RuntimeError(f"BNGsim checkout not found: {bngsim}")

    try:
        inside_work_tree = git(bngsim, ["rev-parse", "--is-inside-work-tree"])
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"not a Git worktree: {bngsim}") from exc
    if inside_work_tree != "true":
        raise RuntimeError(f"not a Git worktree: {bngsim}")

    return bngsim


def verify_clean_source_checkout(bngsim: Path) -> dict[str, Any]:
    changes = git_lines(bngsim, ["status", "--short"])
    if changes:
        preview = "\n  ".join(changes[:10])
        extra = "" if len(changes) <= 10 else f"\n  ... {len(changes) - 10} more"
        raise RuntimeError(
            "BNGsim checkout has local changes. Vendor only from a clean source checkout.\n"
            f"  {preview}{extra}"
        )

    branch = git(bngsim, ["rev-parse", "--abbrev-ref", "HEAD"])
    root = git(bngsim, ["rev-parse", "--show-toplevel"])
    return {
        "repo_root": root,
        "current_branch": branch,
        "origin_remote": remote_url(bngsim, "origin"),
        "head_commit": git(bngsim, ["rev-parse", "HEAD"]),
    }


def resolve_ref(bngsim: Path, ref: str) -> str:
    try:
        return git(bngsim, ["rev-parse", "--verify", f"{ref}^{{commit}}"])
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"could not resolve BNGsim ref '{ref}'") from exc


def file_record(rel_path: str, source_path: str, data: bytes, description: str) -> dict[str, Any]:
    return {
        "path": rel_path,
        "source_path": source_path,
        "sha256": sha256_hex(data),
        "bytes": len(data),
        "description": description,
    }


def bngsim_exprtk_metadata(bngsim: Path, commit: str) -> dict[str, Any] | None:
    try:
        data = git_show(bngsim, commit, "third_party/exprtk/VENDOR.json")
    except subprocess.CalledProcessError:
        return None
    return json_loads_bytes(data)


def build_metadata(
    bngsim: Path,
    checkout_info: dict[str, Any],
    ref: str,
    commit: str,
    source_bytes: dict[str, bytes],
) -> dict[str, Any]:
    files = {}
    for key, spec in FILES.items():
        files[key] = file_record(
            spec["path"], spec["source_path"], source_bytes[key], spec["description"]
        )

    return {
        "name": "BNGsim expression evaluator for RuleMonkey standalone builds",
        "vendored_path": "third_party/bngsim_expr",
        "source": {
            "name": "BNGsim",
            "repo_remote": checkout_info["origin_remote"],
            "repo_subdir": "bngsim",
            "branch_or_ref": ref,
            "current_branch": checkout_info["current_branch"],
            "commit": commit,
        },
        "imported_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "files": files,
        "upstream_exprtk": bngsim_exprtk_metadata(bngsim, commit),
        "local_carries": [],
        "notes": [
            "BNGsim owns the ExprTkEvaluator compatibility layer.",
            "RuleMonkey vendors this layer verbatim for standalone builds.",
            "Inside a BNGsim build, RuleMonkey links the host bngsim::expression target instead.",
            "Refresh these files only through scripts/vendor_exprtk.py.",
        ],
    }


def source_bytes_at_commit(bngsim: Path, commit: str) -> dict[str, bytes]:
    return {
        key: git_show(bngsim, commit, spec["source_path"])
        for key, spec in FILES.items()
    }


def write_vendor_note(metadata: dict[str, Any]) -> None:
    source = metadata["source"]
    lines = [
        "Vendored BNGsim expression evaluator",
        "======================================",
        "",
        "RuleMonkey's standalone rate-law / function / parameter expression",
        "evaluator is bngsim::ExprTkEvaluator. These files are copied verbatim",
        "from BNGsim and are compiled only in standalone RuleMonkey builds;",
        "inside a BNGsim build, RuleMonkey links the host bngsim::expression",
        "target instead.",
        "",
        f"Source repo  : {source.get('repo_remote') or 'unknown'}",
        f"Source subdir: {source.get('repo_subdir')}",
        f"Pinned commit: {source.get('commit')}",
        "Metadata     : third_party/bngsim_expr/VENDOR.json",
        "",
        "Vendored files (copied verbatim, no edits):",
    ]
    for record in metadata["files"].values():
        lines.append(f"  {record['path']} <- bngsim {record['source_path']}")
    lines.extend(
        [
            "",
            "Refresh / drift-check",
            "---------------------",
            "This file is a human-readable companion to VENDOR.json. The JSON",
            "metadata is the source of truth for checksums and provenance.",
            "",
            "  scripts/vendor_exprtk.py --bngsim-repo /path/to/bngsim",
            "  scripts/vendor_exprtk.py --bngsim-repo /path/to/bngsim --summary",
            "  scripts/vendor_exprtk.py --check",
            "",
            "exprtk.hpp is header-only, but version skew between RuleMonkey's",
            "standalone copy and BNGsim's expression layer is undefined behavior",
            "(an ODR violation). The vendoring check guards against that drift.",
        ]
    )
    VENDOR_NOTE.write_text("\n".join(lines) + "\n")


def write_refresh(source_bytes: dict[str, bytes], metadata: dict[str, Any]) -> None:
    for key, spec in FILES.items():
        out = RM_ROOT / spec["path"]
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(source_bytes[key])
    VENDOR_DIR.mkdir(parents=True, exist_ok=True)
    VENDOR_JSON.write_text(json.dumps(metadata, indent=2) + "\n")
    write_vendor_note(metadata)


def metadata_commit(metadata: dict[str, Any] | None) -> str | None:
    if metadata:
        commit = metadata.get("source", {}).get("commit")
        if isinstance(commit, str) and commit:
            return commit
    return read_legacy_pin()


def local_metadata_drift(metadata: dict[str, Any] | None) -> list[str]:
    if metadata is None:
        return [f"missing {VENDOR_JSON}"]

    drift: list[str] = []
    if metadata.get("local_carries") != []:
        drift.append(f"local_carries: expected [], found {metadata.get('local_carries')!r}")

    files = metadata.get("files", {})
    for key, spec in FILES.items():
        record = files.get(key)
        if not isinstance(record, dict):
            drift.append(f"files.{key}: missing metadata")
            continue

        rel_path = record.get("path")
        source_path = record.get("source_path")
        if rel_path != spec["path"]:
            drift.append(f"files.{key}.path: expected {spec['path']!r}, found {rel_path!r}")
            continue
        if source_path != spec["source_path"]:
            drift.append(
                f"files.{key}.source_path: expected {spec['source_path']!r}, "
                f"found {source_path!r}"
            )

        path = RM_ROOT / spec["path"]
        if not path.exists():
            drift.append(f"missing {spec['path']}")
            continue
        data = path.read_bytes()
        actual_sha = sha256_hex(data)
        actual_bytes = len(data)
        if record.get("sha256") != actual_sha:
            drift.append(
                f"files.{key}.sha256: expected {record.get('sha256')!r}, "
                f"found {actual_sha!r}"
            )
        if record.get("bytes") != actual_bytes:
            drift.append(
                f"files.{key}.bytes: expected {record.get('bytes')!r}, "
                f"found {actual_bytes!r}"
            )

    if not metadata_commit(metadata):
        drift.append("source.commit: missing pinned BNGsim commit")
    return drift


def source_drift(bngsim: Path, commit: str) -> list[str]:
    drift: list[str] = []
    for spec in FILES.values():
        source_data = git_show(bngsim, commit, spec["source_path"])
        local_data = (RM_ROOT / spec["path"]).read_bytes()
        if source_data != local_data:
            drift.append(
                f"{spec['path']} differs from BNGsim {commit}:./{spec['source_path']}"
            )
    return drift


def print_summary(bngsim: Path, checkout_info: dict[str, Any], ref: str, commit: str) -> None:
    metadata = read_vendor_metadata()
    current_commit = metadata_commit(metadata) or "missing"
    source_bytes = source_bytes_at_commit(bngsim, commit)
    exprtk_meta = bngsim_exprtk_metadata(bngsim, commit)

    print(f"BNGsim checkout: {bngsim}")
    print(f"BNGsim repo root: {checkout_info['repo_root']}")
    print(f"BNGsim origin: {checkout_info['origin_remote'] or 'unknown'}")
    print(f"BNGsim branch: {checkout_info['current_branch']}")
    print(f"Resolved ref: {ref}")
    print(f"Resolved commit: {commit}")
    print(f"Current RuleMonkey pin: {current_commit}")
    if exprtk_meta:
        header = exprtk_meta.get("files", {}).get("header", {})
        upstream = exprtk_meta.get("source", {})
        print(f"BNGsim upstream ExprTk commit: {upstream.get('commit', 'unknown')}")
        print(f"BNGsim upstream ExprTk sha256: {header.get('sha256', 'unknown')}")

    print("Vendored file deltas:")
    for key, spec in FILES.items():
        local_path = RM_ROOT / spec["path"]
        local_sha = sha256_hex(local_path.read_bytes()) if local_path.exists() else "missing"
        candidate_sha = sha256_hex(source_bytes[key])
        status = "same" if local_sha == candidate_sha else "changed"
        print(f"  {spec['path']}: {status} ({local_sha} -> {candidate_sha})")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--bngsim-repo",
        type=Path,
        help="BNGsim checkout to vendor from; can also be set with BNGSIM_REPO",
    )
    ap.add_argument(
        "--ref",
        default=None,
        help=f"BNGsim ref/tag/commit to vendor on refresh or summarize (default: {DEFAULT_REF})",
    )
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument(
        "--check",
        action="store_true",
        help="verify vendored files against metadata, and against BNGsim if --bngsim-repo is set",
    )
    mode.add_argument(
        "--summary",
        action="store_true",
        help="preview vendoring details for the resolved BNGsim ref; write nothing",
    )
    args = ap.parse_args()

    needs_source = not args.check or args.summary or args.bngsim_repo is not None
    bngsim = resolve_bngsim_repo(args.bngsim_repo, required=needs_source)

    if args.check:
        metadata = read_vendor_metadata()
        drift = local_metadata_drift(metadata)
        if bngsim:
            checkout_info = verify_clean_source_checkout(bngsim)
            del checkout_info  # Source cleanliness is the important guard here.
            commit = metadata_commit(metadata)
            if commit is None:
                drift.append("cannot source-check without a pinned BNGsim commit")
            else:
                drift.extend(source_drift(bngsim, commit))
        if drift:
            print("vendor_exprtk: DRIFT detected:", file=sys.stderr)
            for item in drift:
                print(f"  {item}", file=sys.stderr)
            print(
                "Run scripts/vendor_exprtk.py --bngsim-repo /path/to/bngsim to refresh.",
                file=sys.stderr,
            )
            return 1
        commit = metadata_commit(metadata)
        source_note = " and pinned BNGsim commit" if bngsim else ""
        print(f"vendor_exprtk: OK - vendored files match VENDOR.json{source_note} ({commit})")
        return 0

    assert bngsim is not None
    ref = args.ref or DEFAULT_REF
    checkout_info = verify_clean_source_checkout(bngsim)
    commit = resolve_ref(bngsim, ref)

    if args.summary:
        print_summary(bngsim, checkout_info, ref, commit)
        return 0

    source_bytes = source_bytes_at_commit(bngsim, commit)
    metadata = build_metadata(bngsim, checkout_info, ref, commit, source_bytes)
    write_refresh(source_bytes, metadata)
    print(
        "vendor_exprtk: refreshed 3 files plus VENDOR.json/VENDOR "
        f"from BNGsim {ref} ({commit})"
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        sys.exit(f"vendor_exprtk: {exc}")

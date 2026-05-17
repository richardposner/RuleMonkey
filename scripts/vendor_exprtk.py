#!/usr/bin/env python3
"""Vendor (and drift-check) the ExprTk expression evaluator from BNGsim.

RuleMonkey's rate-law / function / parameter expression evaluator is
bngsim::ExprTkEvaluator (issue #6).  Three files are vendored verbatim from
BNGsim into third_party/; this script refreshes them from a pinned commit
and, with --check, verifies they have not drifted from that pin.

  scripts/vendor_exprtk.py            refresh the 3 files from BNGsim `main`
  scripts/vendor_exprtk.py --ref TAG  ... from a specific BNGsim ref/commit
  scripts/vendor_exprtk.py --check    verify the 3 files still match the pin

The pin is the `Pinned commit:` line of third_party/bngsim_expr/VENDOR.
Drift between RM's copy and BNGsim's is an ODR violation (undefined
behavior) for any binary that ends up with both — hence --check, which is
cheap enough to run in CI on every push.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

RM_ROOT = Path(__file__).resolve().parents[1]
VENDOR_NOTE = RM_ROOT / "third_party" / "bngsim_expr" / "VENDOR"
DEFAULT_BNGSIM = Path("~/Code/PyBNF-Private/bngsim").expanduser()
PIN_PREFIX = "Pinned commit: "

# vendored path (RuleMonkey-relative) -> source path (BNGsim-checkout-relative)
FILES = {
    "third_party/exprtk/exprtk.hpp": "third_party/exprtk/exprtk.hpp",
    "third_party/bngsim_expr/include/bngsim/expression.hpp": "include/bngsim/expression.hpp",
    "third_party/bngsim_expr/src/expression.cpp": "src/expression.cpp",
}


def git_show(bngsim: Path, ref: str, rel: str) -> bytes:
    """Return the bytes of `rel` at `ref`.  `:./path` resolves relative to the
    -C directory, so this works wherever bngsim/ sits inside its git repo."""
    return subprocess.run(
        ["git", "-C", str(bngsim), "show", f"{ref}:./{rel}"],
        check=True,
        stdout=subprocess.PIPE,
    ).stdout


def read_pin() -> str:
    for line in VENDOR_NOTE.read_text().splitlines():
        if line.startswith(PIN_PREFIX):
            return line[len(PIN_PREFIX) :].strip()
    sys.exit(f"vendor_exprtk: no '{PIN_PREFIX}' line in {VENDOR_NOTE}")


def write_pin(commit: str) -> None:
    lines = VENDOR_NOTE.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith(PIN_PREFIX):
            lines[i] = PIN_PREFIX + commit
            break
    else:
        sys.exit(f"vendor_exprtk: no '{PIN_PREFIX}' line in {VENDOR_NOTE}")
    VENDOR_NOTE.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--bngsim-repo",
        type=Path,
        default=DEFAULT_BNGSIM,
        help="BNGsim checkout to vendor from (default: ~/Code/PyBNF-Private/bngsim)",
    )
    ap.add_argument(
        "--ref", default="main", help="BNGsim ref/tag/commit to vendor on refresh (default: main)"
    )
    ap.add_argument(
        "--check",
        action="store_true",
        help="verify the vendored files match the pin; write nothing",
    )
    args = ap.parse_args()

    bngsim = args.bngsim_repo.expanduser()
    if not bngsim.is_dir():
        # Soft skip: a CI job or fresh clone without a BNGsim checkout should
        # not hard-fail here (mirrors scripts/clang-tidy-staged.sh).
        print(f"vendor_exprtk: BNGsim checkout not found at {bngsim}; skipping.", file=sys.stderr)
        return 0

    if args.check:
        pin = read_pin()
        drift = [
            dst
            for dst, src in FILES.items()
            if git_show(bngsim, pin, src) != (RM_ROOT / dst).read_bytes()
        ]
        if drift:
            print(
                f"vendor_exprtk: DRIFT — vendored files differ from pinned BNGsim {pin}:",
                file=sys.stderr,
            )
            for dst in drift:
                print(f"  {dst}", file=sys.stderr)
            print("Run scripts/vendor_exprtk.py to refresh, or restore the files.", file=sys.stderr)
            return 1
        print(f"vendor_exprtk: OK — 3 files match pinned BNGsim {pin}")
        return 0

    commit = subprocess.run(
        ["git", "-C", str(bngsim), "rev-parse", "--verify", f"{args.ref}^{{commit}}"],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip()
    for dst, src in FILES.items():
        out = RM_ROOT / dst
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(git_show(bngsim, commit, src))
    write_pin(commit)
    print(f"vendor_exprtk: refreshed 3 files + re-pinned to BNGsim {commit}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except subprocess.CalledProcessError as exc:
        sys.exit(f"vendor_exprtk: git command failed: {exc}")

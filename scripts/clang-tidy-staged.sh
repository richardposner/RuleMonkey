#!/usr/bin/env bash
# Run clang-tidy on the C++ files passed as arguments (typically the
# staged set, fed in by pre-commit).
#
# The binary comes from pre-commit, which installs the pinned clang-tidy
# wheel named in .pre-commit-config.yaml into the hook's venv and puts it
# on PATH.  So "not found" is a real fault now, not the ordinary state of
# a machine without LLVM, and this exits non-zero rather than skipping.
# It used to skip silently, which is how the gate came to be off on a
# working machine without anything saying so.
#
# One soft path remains: no compile_commands.json under build/ exits 0
# with a hint, because a fresh clone genuinely has none until the first
# `cmake --preset release` (the pre-push hook) has run.  It says plainly
# that nothing was linted.

set -euo pipefail

# Resolve clang-tidy: PATH first, then common Homebrew prefixes (so a
# pre-commit hook running under a non-login shell still finds it on macs
# where llvm is installed via `brew install llvm` and only added to PATH
# from the user's interactive .zshrc).
CLANG_TIDY="${CLANG_TIDY:-}"
if [[ -z "$CLANG_TIDY" ]]; then
  if command -v clang-tidy >/dev/null 2>&1; then
    CLANG_TIDY="$(command -v clang-tidy)"
  else
    for cand in /opt/homebrew/opt/llvm/bin/clang-tidy /usr/local/opt/llvm/bin/clang-tidy; do
      if [[ -x "$cand" ]]; then
        CLANG_TIDY="$cand"
        break
      fi
    done
  fi
fi
if [[ -z "$CLANG_TIDY" ]]; then
  echo "[clang-tidy-staged] clang-tidy not found -- NOT LINTED." >&2
  echo "[clang-tidy-staged] pre-commit is supposed to supply it: check the clang-tidy hook's" >&2
  echo "[clang-tidy-staged] additional_dependencies pin in .pre-commit-config.yaml, and try" >&2
  echo "[clang-tidy-staged] 'pre-commit clean && pre-commit install --install-hooks'." >&2
  echo "[clang-tidy-staged] To run this script outside pre-commit, put a clang-tidy on PATH or" >&2
  echo "[clang-tidy-staged] set CLANG_TIDY=/abs/path/to/clang-tidy." >&2
  exit 1
fi

# macOS SDK 26's libc++ headers use __builtin_clzg / __builtin_ctzg, which
# clang front-ends older than 19 do not have.  Such a clang-tidy does not
# fail to run; it emits thousands of clang-diagnostic-errors from inside
# <algorithm> and <charconv> and never reaches any RM code, which reads
# like the tree is broken.  Name it instead.
ct_major="$("$CLANG_TIDY" --version 2>/dev/null | sed -n 's/.*LLVM version \([0-9][0-9]*\).*/\1/p' | head -1)"
if [[ "$(uname -s)" == "Darwin" && -n "$ct_major" && "$ct_major" -lt 19 ]]; then
  echo "[clang-tidy-staged] clang-tidy $ct_major is too old for this macOS SDK -- NOT LINTED." >&2
  echo "[clang-tidy-staged] SDK 26's libc++ needs a clang front-end >= 19; an older one reports" >&2
  echo "[clang-tidy-staged] the standard library as broken instead of linting your code." >&2
  echo "[clang-tidy-staged] found: $CLANG_TIDY" >&2
  exit 1
fi

# Prefer release (matches the pre-push cmake-configure preset), then fall
# back to any other build dir the developer happens to have around.
build_dir=""
for cand in build/release build/debug build/asan build/profile build/devprof; do
  if [[ -f "$cand/compile_commands.json" ]]; then
    build_dir="$cand"
    break
  fi
done
if [[ -z "$build_dir" ]]; then
  for f in build/*/compile_commands.json; do
    [[ -f "$f" ]] || continue
    build_dir="$(dirname "$f")"
    break
  done
fi
if [[ -z "$build_dir" ]]; then
  echo "[clang-tidy-staged] no compile_commands.json under build/ -- NOT LINTED." >&2
  echo "[clang-tidy-staged] run 'cmake --preset release' once to populate it; until then this" >&2
  echo "[clang-tidy-staged] hook passes without checking anything." >&2
  exit 0
fi

# On macOS, the compile DB is typically produced by Apple Clang against
# the Command Line Tools SDK.  Homebrew's clang-tidy ships its own driver
# and won't infer Apple's sysroot, so it can't find <memory>, <cstdint>,
# etc.  Inject -isysroot explicitly when we're on macOS and xcrun knows
# where the SDK lives.
extra_args=()
if [[ "$(uname -s)" == "Darwin" ]] && command -v xcrun >/dev/null 2>&1; then
  sdk_path="$(xcrun --show-sdk-path 2>/dev/null || true)"
  if [[ -n "$sdk_path" && -d "$sdk_path" ]]; then
    extra_args+=(--extra-arg-before="-isysroot" --extra-arg-before="$sdk_path")
  fi
fi

status=0
for f in "$@"; do
  [[ -f "$f" ]] || continue
  if ! "$CLANG_TIDY" -p "$build_dir" --quiet "${extra_args[@]}" "$f"; then
    status=1
  fi
done
exit "$status"

#!/usr/bin/env bash
# Run clang-tidy on the C++ files passed as arguments (typically the
# staged set, fed in by pre-commit).  Degrades gracefully:
#
#   - if clang-tidy isn't installed, prints a hint and exits 0
#   - if no compile_commands.json is present under build/, prints a hint
#     and exits 0 (the user just needs to run `cmake --preset release`
#     once to populate it)
#
# That intentional softness keeps the hook from blocking commits on a
# fresh clone or on a machine without LLVM installed.  CI is the place to
# turn this into a hard gate.

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
  echo "[clang-tidy-staged] clang-tidy not found; skipping." >&2
  echo "[clang-tidy-staged] install: 'brew install llvm' (macOS) or your distro's clang-tidy package." >&2
  echo "[clang-tidy-staged] or set CLANG_TIDY=/abs/path/to/clang-tidy in your environment." >&2
  exit 0
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
  echo "[clang-tidy-staged] no compile_commands.json under build/; skipping." >&2
  echo "[clang-tidy-staged] run 'cmake --preset release' once to populate it." >&2
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

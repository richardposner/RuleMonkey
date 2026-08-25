# Releasing RuleMonkey

How a RuleMonkey release is cut. The procedure was practice before it was
written down; this file records what the `chore(release):` commits have been
doing since 3.1.2, with the reasoning that decided each step.

Downstream vendoring into BNGsim is a separate, later job and is documented in
the *other* repo, at `bngsim/scripts/RULEMONKEY_VENDORING.md`. Step 7 below is
only the handoff.

## Decide The Version

RuleMonkey is 3.x and versions by what an **embedder or a model can observe**,
not by which Keep a Changelog headings the section happens to carry. The
headings are a bad proxy and the history says so: 3.6.1 is a patch with both a
`### Changed` and an `### Added` block, because both were vendoring and CI
infrastructure that no embedder compiles against or runs.

Ask, in order:

| | Test | Precedent |
|---|---|---|
| **MAJOR** | A public header changes signature — something an embedder compiles against. | None yet in 3.x. |
| **MINOR** | Anything an embedder or a model can observe moves: a new Tier-0 refusal, a changed vector an embedder indexes positionally, a new engine capability. | 3.8.0 (new Tier-0 refusals), 3.10.0 (Added block + a `function_names()` vector that got shorter). |
| **PATCH** | Nothing observable moves. Fixes, performance, build, harness, CI, vendoring. | 3.10.1 (two performance fixes), 3.8.1 (ASan instrumentation + a harness script), 3.6.1 (a vendoring re-pin + a CI drift-guard). |

A new Tier-0 refusal is a **minor**, not a major — 3.8.0 set that precedent and
3.10.0 followed it. A model that loaded under the previous version and now stops
at load is a real break for whoever runs it, but the C++ surface is unchanged,
so an embedder still compiles and links untouched.

The evidence for a patch is worth gathering explicitly, because it is cheap and
it is what the commit message has to argue:

```bash
git diff --stat vX.Y.Z..HEAD -- include/
```

Empty output means no embedder-visible signature moved. Then confirm the
changelog entries themselves claim trajectories are unchanged — the corpus
sweeps say "byte-identical output" when they were run — and that no refusal,
warning or diagnostic was added or removed.

## Version Anchors

Four files carry the version and **there is no single source of truth** — unlike
BNGsim, where `pyproject.toml` derives every other anchor. They must move
together in one commit, or the generated CMake package version disagrees with
the citation metadata:

| File | What moves |
|---|---|
| `CMakeLists.txt` | `project(RuleMonkey VERSION X.Y.Z LANGUAGES CXX)` |
| `pyproject.toml` | `version = "X.Y.Z"` — the `rulemonkey-harness` version, kept in lockstep with the engine |
| `CITATION.cff` | `version:` **and** `date-released:` |
| `CHANGELOG.md` | the section heading, plus two link definitions at the bottom |

`build/<preset>/RuleMonkeyConfigVersion.cmake` is generated from
`CMakeLists.txt` at configure time. Never edit it; check it in step 3.

## 1. Close The Changelog

Turn the open `## [Unreleased]` heading into the release section, leaving an
empty `## [Unreleased]` above it for the next cycle:

```markdown
## [Unreleased]

## [X.Y.Z] — YYYY-MM-DD

### Fixed
```

Then update the two link definitions at the bottom of the file — the compare
link moves to the new tag, and the new tag gets its own release link:

```markdown
[Unreleased]: https://github.com/richardposner/RuleMonkey/compare/vX.Y.Z...HEAD
[X.Y.Z]: https://github.com/richardposner/RuleMonkey/releases/tag/vX.Y.Z
```

Not everything that landed since the last tag needs an entry. A refresh of the
vendored `third_party/bngsim_expr` layer is deliberately **not** changelogged —
inside a BNGsim build CMake links the host `bngsim::expression` target and
RuleMonkey's standalone copy is never compiled, so refreshing its contents is
not a change to RuleMonkey's own behavior. `#58` and `#78` both went unlisted
on that reasoning. A re-pin that changes *policy* rather than contents is
different, and 3.6.1 did list one.

If the section's structure has drifted while it was open — duplicate `### Added`
headings, a block sitting after `### Fixed` — repair it now, in this commit, and
say so in the message. 3.10.0 did exactly that.

## 2. Move The Version Anchors

All four files from the table above, in the same commit as step 1.

## 3. Verify Locally

Four gates. Run them all before opening the PR.

```bash
cmake --preset release && cmake --build --preset release && ctest --preset release --output-on-failure
```

Record the pass count in the commit message — it was 47/47 at 3.10.0 and 49/49
at 3.10.1, and a reader should be able to see the two new cases arrive with the
fixes that brought them.

```bash
grep -m1 'set(PACKAGE_VERSION ' build/release/RuleMonkeyConfigVersion.cmake
```

Must report the new version. This is the check that catches a missed
`CMakeLists.txt` anchor.

**Run it after the `cmake --preset release` above, on the release commit.** That
file is generated only when CMake configures, so it reports whatever the tree
said at the last configure — check it against a build directory last configured
on some other branch and it will happily report the *old* version, or a stale
new one. It is evidence only when the configure that produced it is the one you
just ran here.

```bash
uvx --from "mkdocs-material==9.5.*" mkdocs build --strict
```

**`mkdocs` is deliberately not in `.venv`** — `.github/workflows/docs.yml`
installs it ad hoc, so `mkdocs` is not on `PATH` after activating the project
venv and `uvx` is the way to run it at CI's pin. The build writes `site/`, which
is gitignored.

```bash
.venv/bin/pre-commit run --files CHANGELOG.md CMakeLists.txt pyproject.toml CITATION.cff
```

`pre-commit` *is* in `.venv`. Most hooks skip on this file set — that is
expected, not a misconfiguration.

## 4. Commit And Pull Request

Branch `release/X.Y.Z`, commit subject `chore(release): X.Y.Z`, open a PR
against `main`. Every release since 3.8.0 went through a PR rather than a direct
push.

The commit message is where the release argues for itself, and these messages
are long on purpose — the changelog says what changed, the commit says why the
number moved the way it did. Cover:

- What the section contains, and the version decision **with its evidence** —
  the `include/` diff, the trajectory claims, the absence of new refusals.
- Why not the next tier up, naming the precedent it follows.
- Anything deliberately excluded, and why. Both 3.10.0 and 3.10.1 close by
  noting PyPI is out of scope.
- A `Verified:` line listing all four gates from step 3 with their results.

CI runs nine checks: `build` on macOS/Ubuntu/Windows, `asan` on macOS/Ubuntu,
`clang_tidy`, `feature_coverage`, `perf_diff`, and `vendor_check`. Wait for all
nine. `vendor_check` clones `lanl/bngsim` and byte-compares the vendored
expression layer against it, so a green one here is also the first half of the
evidence that step 7 will not be blocked.

Squash-merge.

## 5. Tag

Annotated tag, on the squashed release commit on `main`:

```bash
git checkout main && git pull --ff-only origin main
git tag -a vX.Y.Z -m "RuleMonkey X.Y.Z" <release-commit>
git push origin vX.Y.Z
```

Every 3.x tag is annotated with the message `RuleMonkey X.Y.Z`. A lightweight
tag would break that consistency and lose the tagger date the release list
sorts on.

## 6. Publish The Release

The release notes are the changelog section **verbatim** — not a rewrite, not a
summary:

```bash
awk '/^## \[X.Y.Z\]/{f=1} /^## \[<previous>\]/{f=0} f' CHANGELOG.md > /tmp/relnotes.md
gh release create vX.Y.Z --title "RuleMonkey X.Y.Z" --notes-file /tmp/relnotes.md
```

No assets. RuleMonkey is consumed as source — vendored into BNGsim, or built
from a checkout — so there is nothing to attach and no release has ever carried
one.

Confirm it landed as `Latest`:

```bash
gh release list --limit 3
```

## 7. Downstream: Re-Vendor Into BNGsim

A RuleMonkey release does not reach BNGsim on its own. That refresh is run from
the BNGsim repo and is documented there, in
`bngsim/scripts/RULEMONKEY_VENDORING.md`. Two things about it are worth knowing
from this side:

- **It may need work landed here first.** The `EXPRTK_SYNC_FILES` guard compares
  RuleMonkey's standalone `third_party/bngsim_expr` copy against BNGsim's tree
  and fails closed. If BNGsim's expression layer has moved since that copy was
  last synced, `scripts/vendor_exprtk.py --bngsim-repo <bngsim>` has to run here
  and land upstream *before* the vendor refresh can run at all. That is what
  `#78` was for, and it is why `#58` exists as a standalone refresh commit rather
  than a release.
- **No release tag is required.** BNGsim vendors `main` at a commit, not a tag.
  The pin before 3.10.0 was `1d14160`, an ordinary refresh commit.

Whether the ordered two-repo dance is needed is visible before any writing
happens, from the `ExprTk sync pin` line of:

```bash
python3 bngsim/scripts/vendor_rulemonkey.py --rulemonkey-repo /tmp/rulemonkey-vendor-candidate --ref main --summary
```

## Not In Scope: PyPI

Neither `rulemonkey` nor `rulemonkey-harness` exists on PyPI, and
`pyproject.toml` is not a distribution — it declares no build backend, and
`harness/` is a flat directory of scripts with no package layout. There is
nothing to upload. Tracked separately; both 3.10.0 and 3.10.1 say so explicitly
in their release commits.

#!/usr/bin/env bash
# =============================================================================
# scripts/release.sh  --  release a new code bank version (and hand it to the website)
# =============================================================================
# Usage, from the code bank folder, after committing your changes:
#   scripts/release.sh 0.6.1                 check, test, tag v0.6.1, start the website's release check
#   scripts/release.sh 0.6.1 --quick         skip the long pipeline tests (only the manifest and changed modules)
#   scripts/release.sh 0.6.1 --no-site      only tag; do not start the website's release check
#
# Before tagging it checks that:
#   - everything is committed
#   - CODE_BANK_VERSION in R/utils.R, `version:` in manifest/tools.yaml and a "## v0.6.1" CHANGELOG entry agree
#   - tests/check_manifest.R passes (the manifest matches the defaults and only uses known site building blocks)
#   - the tests of every module changed since the last release pass (and the main pipeline tests, unless --quick)
# The tag records who released it ("released-by: $USER"). That person approves the website change:
# the website ($COPILOT_HOME/copilot; COPILOT_HOME defaults to the folder above this code bank) checks itself
# against the release in a staging copy, lets its release agent fix what is needed, and writes a review; then
#   $COPILOT_HOME/copilot/ops/release approve v0.6.1
# The check runs containers/ada/release.slurm (with the coding model): as a SLURM job where there is SLURM
# (RELEASE_WITH_SLURM=0: here instead), otherwise here in the foreground (fuji).
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
VERSION="${1:?usage: scripts/release.sh <version, e.g. 0.6.1> [--quick] [--no-site]}"
VERSION="${VERSION#v}"
TAG="v$VERSION"
QUICK=0; SITE=1
for a in "${@:2}"; do
  case "$a" in --quick) QUICK=1 ;; --no-site) SITE=0 ;; *) echo "Unknown option $a"; exit 1 ;; esac
done
RSCRIPT="${RSCRIPT:-Rscript}"
fail() { echo "RELEASE STOPPED: $*" >&2; exit 1; }

[ -z "$(git status --porcelain --untracked-files=no)" ] || fail "commit your changes first (git status shows uncommitted changes)."
git rev-parse -q --verify "refs/tags/$TAG" >/dev/null && fail "$TAG already exists. Use a new version number."
grep -q "CODE_BANK_VERSION <- \"$VERSION\"" R/utils.R || fail "R/utils.R must say CODE_BANK_VERSION <- \"$VERSION\"."
grep -q "^version: $VERSION\b" manifest/tools.yaml || fail "manifest/tools.yaml must say version: $VERSION."
grep -q "^## $TAG\b" CHANGELOG.md || fail "CHANGELOG.md needs a '## $TAG' section saying what changed."

echo "== Manifest"
"$RSCRIPT" tests/check_manifest.R

PREV="$(git describe --tags --abbrev=0 2>/dev/null || true)"
CHANGED="$( [ -n "$PREV" ] && git diff --name-only "$PREV"..HEAD || git ls-files )"
for mod in splicing isoforms compare later; do
  if echo "$CHANGED" | grep -q "^$mod/\|^R/"; then
    t=$(ls "$mod"/tests/run_*tests.R 2>/dev/null | head -1)
    if [ -n "$t" ]; then
      if [ "$mod" = later ] && ! command -v STAR >/dev/null; then echo "== $mod tests skipped (STAR, salmon, rMATS not installed here)"; continue; fi
      echo "== $mod tests ($t)"; "$RSCRIPT" "$t"
    fi
  fi
done
if [ "$QUICK" = 0 ] && echo "$CHANGED" | grep -q "^scripts/\|^R/\|^config/"; then
  echo "== Main pipeline tests"; "$RSCRIPT" tests/run_tests.R
fi

git tag -a "$TAG" -m "Code bank $TAG
released-by: ${USER:-$(id -un)}"
echo "Tagged $TAG (released by ${USER:-$(id -un)})."

if [ "$SITE" = 1 ]; then
  COPILOT_HOME="${COPILOT_HOME:-$(cd .. && pwd)}"
  SITE_DIR="$COPILOT_HOME/copilot"
  if [ -x "$SITE_DIR/ops/release" ]; then
    echo "== Website release check ($SITE_DIR)"
    if [ "${RELEASE_WITH_SLURM:-1}" = 1 ] && command -v sbatch >/dev/null && [ -f "$SITE_DIR/containers/ada/release.slurm" ]; then
      (cd "$COPILOT_HOME" && sbatch --export=ALL,TAG="$TAG",COPILOT_HOME="$COPILOT_HOME",COPILOT_CODEBANK_REPO="$(pwd)" \
        "$SITE_DIR/containers/ada/release.slurm")
      echo "Submitted (log: $COPILOT_HOME/rnaseq-release-<job>.log). When it finishes, read"
      echo "$COPILOT_HOME/staging/$TAG/RELEASE_REVIEW.md and approve with: $SITE_DIR/ops/release approve $TAG"
    elif [ -f "$SITE_DIR/containers/ada/release.slurm" ] && [ -f "$COPILOT_HOME/rnaseq-copilot.sif" ]; then
      # no SLURM (e.g. fuji): the same job, here and now (starts or reuses the model server on the GPU)
      TAG="$TAG" COPILOT_HOME="$COPILOT_HOME" COPILOT_CODEBANK_REPO="$(pwd)" bash "$SITE_DIR/containers/ada/release.slurm" || true
    else
      COPILOT_HOME="$COPILOT_HOME" COPILOT_CODEBANK_REPO="$(pwd)" "$SITE_DIR/ops/release" integrate "$TAG" || true
    fi
  else
    echo "No website found at $SITE_DIR (set COPILOT_HOME); the release is tagged only."
  fi
fi

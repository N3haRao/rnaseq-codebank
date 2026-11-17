#!/usr/bin/env bash
# =============================================================================
# scripts/rscript_image.sh  --  Rscript from the site image (rnaseq-copilot.sif)
# =============================================================================
# For servers whose own R lacks the code bank's packages (e.g. fuji): runs the
# image's Rscript with $COPILOT_HOME bound at the same path and the current folder
# kept, so `scripts/release.sh` and the module tests behave as on a full install.
# scripts/release.sh uses it by itself when the server's Rscript is missing packages.
#   COPILOT_HOME=~/rnaseq-copilot scripts/rscript_image.sh tests/check_manifest.R
# =============================================================================
set -euo pipefail
COPILOT_HOME="${COPILOT_HOME:-$(cd "$(dirname "$0")/../.." && pwd)}"
SIF="${COPILOT_SIF:-$COPILOT_HOME/rnaseq-copilot.sif}"
A=$(command -v apptainer || command -v singularity || true)
[ -n "$A" ] || { echo "rscript_image.sh: no apptainer or singularity on this server." >&2; exit 127; }
[ -f "$SIF" ] || { echo "rscript_image.sh: no site image at $SIF." >&2; exit 127; }
mkdir -p "$COPILOT_HOME/tmp"
BINDS="$COPILOT_HOME,$COPILOT_HOME/tmp:/tmp"
case "$PWD" in "$COPILOT_HOME"|"$COPILOT_HOME"/*) ;; *) BINDS="$BINDS,$PWD" ;; esac
exec "$A" exec --containall --bind "$BINDS" --env "COPILOT_HOME=$COPILOT_HOME" --pwd "$PWD" \
  "$SIF" /opt/conda/envs/r/bin/Rscript "$@"

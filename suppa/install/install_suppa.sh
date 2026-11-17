#!/usr/bin/env bash
# =============================================================================
# suppa/install/install_suppa.sh  --  install SUPPA2 for the suppa/ module (once per server)
# =============================================================================
# Makes a conda environment with SUPPA2 at $COPILOT_HOME/envs/suppa (or the folder
# given as $1). The module finds it there by itself; elsewhere set suppa_command
# (or SUPPA_COMMAND) to <env>/bin/suppa.py.
#
# SUPPA2 is MIT-licensed. It needs older Python libraries: with pandas >= 2 or
# statsmodels >= 0.14 diffSplice fails (statsmodels moved multipletests), so the
# versions are pinned below.
#
# Uses micromamba: the one inside the site image (rnaseq-copilot.sif) if there is
# one, else micromamba / mamba / conda on the PATH.
#
#   suppa/install/install_suppa.sh                    # -> $COPILOT_HOME/envs/suppa
#   suppa/install/install_suppa.sh /path/to/envs/suppa
# =============================================================================
set -euo pipefail
COPILOT_HOME=${COPILOT_HOME:-$HOME/rnaseq-copilot}
ENV=${1:-$COPILOT_HOME/envs/suppa}
SPECS=(python=3.10 "suppa>=2.3,<2.5" "pandas<2" "numpy<1.25" "scipy<1.11" "statsmodels=0.13.5" "scikit-learn<1.3")

if [ -x "$ENV/bin/suppa.py" ]; then echo "SUPPA2 is already installed in $ENV"; exit 0; fi
mkdir -p "$(dirname "$ENV")"
if [ -f "$COPILOT_HOME/rnaseq-copilot.sif" ] && { command -v apptainer >/dev/null || command -v singularity >/dev/null; }; then
  A=$(command -v apptainer || command -v singularity)
  echo "Installing SUPPA2 into $ENV with the site image's micromamba ..."
  "$A" exec --bind "$(dirname "$ENV")" "$COPILOT_HOME/rnaseq-copilot.sif" micromamba create -y \
    -r "$(dirname "$ENV")/.mamba" -p "$ENV" -c conda-forge -c bioconda "${SPECS[@]}"
elif command -v micromamba >/dev/null; then
  micromamba create -y -p "$ENV" -c conda-forge -c bioconda "${SPECS[@]}"
elif command -v mamba >/dev/null; then
  mamba create -y -p "$ENV" -c conda-forge -c bioconda "${SPECS[@]}"
elif command -v conda >/dev/null; then
  conda create -y -p "$ENV" -c conda-forge -c bioconda "${SPECS[@]}"
else
  echo "No micromamba, mamba or conda found (and no site image)." >&2; exit 1
fi
rm -rf "$(dirname "$ENV")/.mamba/pkgs"
"$ENV/bin/python" "$ENV/bin/suppa.py" -h | head -3
echo "SUPPA2 installed: $ENV/bin/suppa.py"
echo "Test the module:  SUPPA_ENV=$ENV Rscript suppa/tests/run_suppa_tests.R   (from the code bank folder)"

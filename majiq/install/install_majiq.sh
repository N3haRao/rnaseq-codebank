#!/usr/bin/env bash
# =============================================================================
# majiq/install/install_majiq.sh  --  install MAJIQ v3 for the majiq/ module (once per server)
# =============================================================================
# MAJIQ is NOT open source: it is free for non-commercial use by individuals at academic or
# not-for-profit institutions under its academic licence. Register and accept the licence first:
#   https://majiq.biociphers.org/app_download/
# This script is therefore not run by the site image; each server installs MAJIQ under the lab's
# licence. The module finds it at $COPILOT_HOME/envs/majiq/bin and the licence file at
# $COPILOT_HOME/majiq_license*.lic (or set majiq_bin_dir / majiq_license_file).
#
# What it does (the steps the lab used for MAJIQ 3.0.23, made non-interactive):
#   1. the academic licence file -> $COPILOT_HOME/majiq_license_academic_official.lic
#   2. a conda environment (Python 3.12, htslib, compilers, compression libraries, cmake, pybind11)
#   3. pip install MAJIQ v3 + voila from the academic repository (compiles C++ code, ~10-20 min)
#
#   majiq/install/install_majiq.sh                 # -> $COPILOT_HOME/envs/majiq
#   majiq/install/install_majiq.sh /path/to/env
# Docs: https://biociphers.bitbucket.io/majiq-docs
# =============================================================================
set -euo pipefail
COPILOT_HOME=${COPILOT_HOME:-$HOME/rnaseq-copilot}
E=${1:-$COPILOT_HOME/envs/majiq}
LIC=$COPILOT_HOME/majiq_license_academic_official.lic
mkdir -p "$(dirname "$E")" "$COPILOT_HOME"

## ---- 1. Licence file (MAJIQ >= 2.5 does not run without it) ----
if ! ls "$COPILOT_HOME"/majiq_license* >/dev/null 2>&1; then
  echo "Downloading the academic licence file (you must have registered and accepted the licence terms) ..."
  curl -sSfL -o "$LIC" https://majiq.biociphers.org/app_download/majiq_license_academic_official.lic
fi

## ---- 2. Conda environment ----
SPECS=(python=3.12 htslib zlib pip git "gcc_linux-64>=9" "gxx_linux-64>=9" bzip2 xz zstd cmake "pybind11=2.13.6")
if [ ! -x "$E/bin/python" ]; then
  if [ -f "$COPILOT_HOME/rnaseq-copilot.sif" ] && { command -v apptainer >/dev/null || command -v singularity >/dev/null; }; then
    A=$(command -v apptainer || command -v singularity)
    "$A" exec --bind "$(dirname "$E")" "$COPILOT_HOME/rnaseq-copilot.sif" micromamba create -y \
      -r "$(dirname "$E")/.mamba" -p "$E" -c conda-forge -c bioconda --override-channels "${SPECS[@]}"
  elif command -v micromamba >/dev/null; then micromamba create -y -p "$E" -c conda-forge -c bioconda --override-channels "${SPECS[@]}"
  elif command -v mamba >/dev/null; then mamba create -y -p "$E" -c conda-forge -c bioconda --override-channels "${SPECS[@]}"
  elif command -v conda >/dev/null; then conda create -y -p "$E" -c conda-forge -c bioconda --override-channels "${SPECS[@]}"
  else echo "No micromamba, mamba or conda found (and no site image)." >&2; exit 1; fi
  rm -rf "$(dirname "$E")/.mamba/pkgs"
fi

## ---- 3. MAJIQ v3 + voila from the academic repository ----
## Only the environment's tools: a system-wide anaconda would otherwise give CMake the wrong
## libraries (BZip2::BZip2 not found); SETUPTOOLS_SCM_PRETEND_VERSION avoids a metadata mismatch
## when building rna_moccasin.
if [ ! -x "$E/bin/majiq" ]; then
  env -u CMAKE_PREFIX_PATH PATH="$E/bin:/usr/bin:/bin" \
    CC="$E/bin/x86_64-conda-linux-gnu-gcc" CXX="$E/bin/x86_64-conda-linux-gnu-g++" \
    HTSLIB_LIBRARY_DIR="$E/lib" HTSLIB_INCLUDE_DIR="$E/include" CMAKE_PREFIX_PATH="$E" \
    CMAKE_ARGS="-DCMAKE_PREFIX_PATH=$E -Dpybind11_DIR=$E/share/cmake/pybind11" \
    LDFLAGS="-L$E/lib -Wl,-rpath,$E/lib" CPPFLAGS="-I$E/include" \
    SETUPTOOLS_SCM_PRETEND_VERSION=3.0.23 \
    "$E/bin/pip" install --no-cache-dir "git+https://bitbucket.org/biociphers/majiq_academic.git"
fi

## ---- 4. Check ----
"$E/bin/majiq" --version
"$E/bin/voila" --help | head -3
[ -x "$E/bin/majiq-v3" ] || echo "NOTE: no majiq-v3 command: set run_modulize: false (voila modulize needs majiq-v3 sg-coverage)."
echo "MAJIQ installed in $E; licence: $(ls "$COPILOT_HOME"/majiq_license* | head -1)"
echo "The GFF3 annotation must be the same release as the GTF used for alignment, e.g. for GENCODE 44:"
echo "  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gff3.gz"

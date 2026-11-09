#!/usr/bin/env bash
# ============================================================
# Install R 4.6 + Bioconductor 3.23 binaries on a Colab VM via r2u.
#
# r2u (eddelbuettel/r2u) ships CRAN and Bioconductor packages as
# Ubuntu .deb binaries and wires install.packages() to apt (bspm), so
# DESeq2/minfi/clusterProfiler install in seconds instead of ~40 min
# of compiling. Its current build is R 4.6 / Bioconductor 3.23 — the
# same versions the v4 results were produced with locally.
#
# Safe to rerun: skips work already done on this VM.
# ============================================================
set -euo pipefail

CODENAME="$(lsb_release -cs)"
echo "Ubuntu codename: ${CODENAME}"

if [ -f /etc/apt/sources.list.d/cranapt.list ] && R --version 2>/dev/null | grep -q "R version 4.6"; then
  echo "r2u + R 4.6 already set up on this VM."
else
  case "${CODENAME}" in
    jammy|noble|resolute) ;;
    *) echo "r2u does not support '${CODENAME}'. Falling back to source installs (slow)."; exit 0 ;;
  esac
  wget -q -O /tmp/add_cranapt.sh \
    "https://raw.githubusercontent.com/eddelbuettel/r2u/master/inst/scripts/add_cranapt_${CODENAME}.sh"
  bash /tmp/add_cranapt.sh
fi

# System libraries needed by any package that still builds from source
apt-get install -y -qq --no-install-recommends \
  libcurl4-openssl-dev libssl-dev libxml2-dev libglpk-dev libgsl-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
  libpng-dev libjpeg-dev libtiff-dev zlib1g-dev libbz2-dev liblzma-dev \
  > /dev/null

R --version | head -1
echo "R setup done."

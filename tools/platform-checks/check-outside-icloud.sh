#!/bin/zsh
# Build and check somewhere the file synchroniser is not watching.
#
# The working copy sits under ~/Desktop, which iCloud Drive syncs, and the
# synchroniser writes conflict copies -- "tilting 2.R" beside "tilting.R" --
# into directories while R is building in them. One such copy has already
# appeared in R/, where it would have shipped in the tarball and been sourced
# at install time, and another appeared inside .Rcheck during a check, which
# R CMD check then reported as a non-standard file. Building elsewhere removes
# the interference rather than working around it.
set -e
PKG=${1:-$(pwd)}
WORK=$(mktemp -d /tmp/drbayes-check.XXXXXX)
trap 'echo "left in $WORK"' EXIT

( cd "$PKG" && R CMD build . )
mv "$PKG"/DRBayes_*.tar.gz "$WORK"/
cd "$WORK"
_R_CHECK_FORCE_SUGGESTS_=false _R_CHECK_LIMIT_CORES_=TRUE \
  R CMD check --as-cran DRBayes_*.tar.gz

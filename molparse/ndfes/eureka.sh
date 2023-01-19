#!/bin/bash

ARGS=$@

module load MPICH2

FEHOME=$HOME/tars_n_zips/fe-toolkit
FEBIN=$FEHOME/local/bin

export PREFIX=$FEHOME/local
export LIBDIR=${PREFIX}/lib
export LIB64DIR=${PREFIX}/lib64
export BINDIR=${PREFIX}/bin
export INCDIR=${PREFIX}/include
export LDFLAGS="-L${LIBDIR} -L${LIB64DIR}"
export LIBRARY_PATH="${LIBDIR}:${LIB64DIR}:${LIBRARY_PATH}"
export LD_LIBRARY_PATH="${LIBDIR}:${LIB64DIR}:${LD_LIBRARY_PATH}"
export CPATH="${NUMPYINC}:${PREFIX}/include:${CPATH}"
export PYTHONPATH="${LIBDIR}/python${PYTHONVER}/site-packages:${PYTHONPATH}"

_CONDA_ROOT="/users/mw00368/tars_n_zips/fe-toolkit/local"
# Copyright (C) 2012 Anaconda, Inc
# SPDX-License-Identifier: BSD-3-Clause
\. "$_CONDA_ROOT/etc/profile.d/conda.sh" || return $?
conda activate

$FEBIN/ndfes $ARGS

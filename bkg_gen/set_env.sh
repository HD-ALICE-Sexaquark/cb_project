#!/bin/bash

export PYTHIA_ROOT=$(readlink -f pythia8309_root)
export PYTHIA8DATA=${PYTHIA_ROOT}/share/Pythia8/xmldoc
export PYTHIA8=${PYTHIA_ROOT}

export PATH=${PYTHIA_ROOT}/bin:${PATH}
export LD_LIBRARY_PATH=${PYTHIA_ROOT}/lib:${LD_LIBRARY_PATH}
export ROOT_INCLUDE_PATH=${PYTHIA_ROOT}/include:${ROOT_INCLUDE_PATH}

echo "Setting environment to use PYTHIA8/ANGANTYR"
echo "==========================================="
echo ">> PYTHIA_ROOT       = ${PYTHIA_ROOT}"
echo ">> PYTHIA8DATA       = ${PYTHIA8DATA}"
echo ">> PYTHIA8           = ${PYTHIA8}"
echo ">> PATH              = ${PATH}"
echo ">> LD_LIBRARY_PATH   = ${LD_LIBRARY_PATH}"
echo ">> ROOT_INCLUDE_PATH = ${ROOT_INCLUDE_PATH}"

#!/bin/bash

BKG_PDG_CODE=${1}

CURRENT_DIR=${PWD}

# define build dir
BUILD_DIR=$(readlink -f ./B2a_CB_build)

cd ${BUILD_DIR}

./exampleB2a ${BKG_PDG_CODE}

cd ${CURRENT_DIR}

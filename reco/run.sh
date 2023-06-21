#!/bin/bash

CURRENT_DIR=${PWD}

# define build dir
BUILD_DIR=$(readlink -f ./B2a_MK_build)

cd ${BUILD_DIR}

./exampleB2a
# ./exampleB2a ../../output/run000/event000_sig.csv ../../output/run000/event000_bkg.csv

cd ${CURRENT_DIR}

#!/bin/bash

CURRENT_DIR=${PWD}

# define build dir
BUILD_DIR=$(readlink -f ./B2a_CB_build)

cd ${BUILD_DIR}

./exampleB2a

cd ${CURRENT_DIR}
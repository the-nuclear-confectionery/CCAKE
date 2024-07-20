#!/bin/bash
#----------

BUILD_DIR=build

# Uncomment to build with g++
mkdir -p ${BUILD_DIR}
( cd ${BUILD_DIR} \
  && cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc ..     \
  && make         \
  && make install )

# Uncomment to build with clang++
# mkdir -p ${BUILD_DIR}
# ( cd ${BUILD_DIR} \
#   && cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang ..     \
#   && make         \
#   && make install )

#!/bin/bash

mkdir -p ../Debug/

cd ../Debug/
cmake \
 -DCMAKE_BUILD_TYPE=Debug \
 ../Code/

mkdir -p ../Release/
cd ../Release/
cmake \
  -DCMAKE_INSTALL_PREFIX=$HOME \
  -DCMAKE_BUILD_TYPE=Release \
  ../Code/

cd ../Code/


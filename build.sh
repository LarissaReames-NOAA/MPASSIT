#! /usr/bin/env bash
#
# Author: Larissa Reames CIWRO/NOAA/NSSL/FRDD  

set -eux

target=${target:-"NULL"}
compiler=${compiler:-"intel"}

if [[ "$target" == "linux.*" || "$target" == "macosx.*" ]]; then
 unset -f module
 set +x
 source ./modulefiles/build.$target > /dev/null
 set -x
else
 set +x
 source ./machine-setup.sh
 module use ./modulefiles
 module load build.$target.$compiler.lua > /dev/null
 module list
 set -x
fi

if [[ "$target" == "hera" || "$target" == "orion" || "$target" == "wcoss2" ]]; then
  #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Debug"
   CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF"
   #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DENABLE_DOCS=ON -DBUILD_TESTING=ON"
else
  CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Debug"
  #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF"
fi

rm -fr ./build
mkdir ./build && cd ./build

cmake .. ${CMAKE_FLAGS}

make -j 8 VERBOSE=1
make install

#make test
#ctest -I 4,5

exit

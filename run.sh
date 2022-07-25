#!/bin/bash

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/gitrepos/gmp-6.2.1:$HOME/gitrepos/mpfr-4.1.0:$HOME/gitrepos/next-posits/fftw-3.3.10
export MPC_DIR=$PWD/mpc-1.2.1
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/Users/sameer/gitrepos/next-posits/build

mkdir build
cd build
cmake ../src -DCMAKE_INSTALL_PREFIX=$PWD -DENABLE_FFTW=ON \
      -DFFTW_INCLUDES=/Users/sameer/gitrepos/tfhe/fftw-3.3.10/include \
      -DFFTW_LIBRARIES=/Users/sameer/gitrepos/tfhe/fftw-3.3.10/lib/libfftw3.a \
      -DCMAKE_C_FLAGS="-I$MPC_DIR/include" \
      -DCMAKE_EXE_LINKER_FLAGS="$MPC_DIR/lib/libmpc.a"

make -j
make install

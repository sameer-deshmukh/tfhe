#!/bin/bash

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HOME/gitrepos/gmp-6.2.1:$HOME/gitrepos/mpfr-4.1.0:$HOME/gitrepos/next-posits/fftw-3.3.10
export MPC_DIR=$HOME/gitrepos/next-posits/mpc-1.2.1
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/Users/sameer/gitrepos/next-posits/build

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/sameer/gitrepos/next-posits/build/lib

rm -rf build
mkdir build
cd build
cmake ../src -DCMAKE_INSTALL_PREFIX=$PWD -DENABLE_FFTW=ON -DENABLE_NAYUKI_AVX=OFF -DENABLE_NAYUKI_PORTABLE=OFF \
      -DENABLE_SPQLIOS_AVX=OFF -DENABLE_SPQLIOS_FMA=OFF -DENABLE_TESTS=ON \
      -DFFTW_INCLUDES=/Users/sameer/gitrepos/tfhe/fftw-3.3.10/include \
      -DFFTW_LIBRARIES=/Users/sameer/gitrepos/tfhe/fftw-3.3.10/lib/libfftw3.a \
      -DCMAKE_CXX_FLAGS="-I$MPC_DIR/include" \
      -DCMAKE_EXE_LINKER_FLAGS="$MPC_DIR/lib/libmpc.a /Users/sameer/gitrepos/next-posits/build/lib/libnext_posits.a"

make VERBOSE=1
make install

#!/bin/bash

mkdir -p build
cd build

gcc_path=$(which gcc)
gxx_path=$(which g++)

if [ -z "$gcc_path" ]; then
    echo "gcc not found in PATH"
else
    echo "gcc found at: $gcc_path"
fi

if [ -z "$gxx_path" ]; then
    echo "g++ not found in PATH"
else
    echo "g++ found at: $gxx_path"
fi

export CC=$gcc_path
export CXX=$gxx_path

CPU_CORES=$(nproc)
echo "Building with $CPU_CORES cores."

cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$CPU_CORES

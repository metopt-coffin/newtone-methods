#!/bin/sh

[ ! -d build ] && mkdir build
cd build

cmake \
    -G "Unix Makefiles" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DCMAKE_CXX_COMPILER=clang++ \
    ..

job_cnt=$(lscpu | grep '^CPU(s):' | sed "s/CPU(s)://g")

make -sj `expr $job_cnt / 2`

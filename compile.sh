#!/bin/bash

build_dir="lib"
src_dir="src"

[[ -e $build_dir ]] || mkdir $build_dir
gcc ${src_dir}/lib/{random.c,simulat.c} \
    -I ${src_dir}/include \
    -lm \
    --shared \
    -fPIC \
    -o ${build_dir}/libsimulat.so \

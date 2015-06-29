#!/bin/bash

build_dir="bin"
src_dir="src"

[[ -e $build_dir ]] || mkdir $build_dir
gcc ${src_dir}/lib/{random.c,simulat.c} \
    -I ${src_dir}/include \
    -lm \
    --shared \
    -fPIC \
    -Wall \
    -o ${build_dir}/libsimulat.so \

gcc ${src_dir}/test/test_BPSK.c ${src_dir}/lib/{simulat.c,random.c} \
    -I ${src_dir}/include \
    -lm \
    -lfftw3 \
    -Wall \
    -o ${build_dir}/bpsk

gcc ${src_dir}/test/test_QPSK.c ${src_dir}/lib/{simulat.c,random.c} \
    -I ${src_dir}/include \
    -lm \
    -lfftw3 \
    -Wall \
    -o ${build_dir}/qpsk

gcc ${src_dir}/test/test_OFDM.c ${src_dir}/lib/{simulat.c,random.c}  \
    -I ${src_dir}/include \
    -lm \
    -lfftw3 \
    -Wall \
    -o ${build_dir}/ofdm

#ifndef SIMULAT_H
#define SIMULAT_H

#include <stdbool.h>
#include <fftw3.h>
#include <math.h>

typedef unsigned long ulong;

typedef struct Point {
    double SNR_db;
    double BER;
} Point;

/*串行序列*/
typedef struct Serial {
    double * signal;
    ulong N;
} Serial;

/*串行序列, but存储的是复数*/
typedef struct Serial_C {
    fftw_complex * signal;
    ulong N;
} Serial_C;

/*m路并行*/
typedef struct Parallel {
    double ** signal;
    ulong M;
    ulong N;
} Parallel;
/*m路并行, but存储的是复数*/
typedef struct Parallel_C {
    fftw_complex ** signal;
    ulong M;
    ulong N;
} Parallel_C;

double gen_uniform_random(void);
double gen_rayleigh_random(double sigma);
double gen_standard_normal_random(void);
double gen_normal_random(double mean, double sigma);
int gen_binomial_random(double p);

Serial* gen_signal(ulong N);
void bipolar(Serial* serial, double scale);
void add_noise(Serial * serial, double mean, double sigma);
void add_noise_c(Serial_C * serial_c, double mean, double sigma);
void add_noise_and_fading(Serial * serial, double n_mean,
        double n_sigma, double r_sigma);
void add_noise_and_fading_c(Serial_C * serial, double n_mean,
        double n_sigma, double r_sigma);
Parallel* serial_to_parallel(Serial * serial, ulong M);//因为要补零,所以没有加const
Parallel_C* serial_to_parallel_c(Serial_C* serial, ulong M);
Serial_C* parallel_to_serial_c(Parallel_C* parallel_c);
void add_cycle_prefix(Parallel_C* parallel_c);
void del_cycle_prefix(Parallel_C* parallel_c);
Serial* BPSK(Serial* serial);
Serial* rBPSK(Serial* serial);
Serial_C* QPSK(Serial* serial);
Serial* rQPSK(Serial_C* serial_c);
Serial_C* QAM16(Serial* serial);
Serial* rQAM16(Serial_C* serial_c);
void judge4QAM16(Serial_C* serial_c);
void judge4QPSK(Serial_C* serial_c);
void judge4BPSK(Serial* serial);
void fft(Parallel_C* parallel_c);
void ifft(Parallel_C* parallel_c);
void transpose_array(Parallel* parallel);//转置一个实数矩阵

#endif

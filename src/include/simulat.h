#ifndef SIMULAT_H
#define SIMULAT_H

#include <stdbool.h>

typedef unsigned long ulong;
typedef unsigned short ushort;

typedef struct Point {
    double SNR_db;
    double BER;
} Point;

/*串行*/
typedef struct Serial {
    double * signal;
    ulong N;
} Serial;

/*m路并行*/
typedef struct Parallel {
    double ** signal;
    ushort M;
    ulong N;
} Parallel;

Serial* gen_signal(ulong N);
void free_serial(Serial * serial);
void add_noise(Serial * serial, double mean, double sigma);
void add_noise_and_fading(Serial * serial, double n_mean,
        double n_sigma, double r_sigma);
Parallel* serial_to_parallel(Serial * serial, ushort M);
void free_parallel(Parallel * parallel);
double gen_uniform_random(void);
double gen_rayleigh_random(double sigma);
double gen_standard_normal_random(void);
double gen_normal_random(double mean, double sigma);
int gen_binomial_random(double p);


/*
 *参数:
 *    N:         产生的信号比特数
 *    DB_MAX:    SNR_db最大值
 *    fading:    是否添加瑞利衰落
 *    SNR_BER_p: 信噪比/误码率--作为返回值, 指向Point结构组成的数组
 */

void simulat_BPSK(const ulong N, const ushort DB_MAX,
        const bool fading, Point* SNR_BER_p);

#endif

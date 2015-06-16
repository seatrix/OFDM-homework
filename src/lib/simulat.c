#include "simulat.h"

#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>


void simulat_BPSK(const ulong N, const ushort DB_MAX,
        const bool fading, Point* snr_ber_p)
{
    /*先计算噪声标准差*/
    double sigma[DB_MAX];
    for (int i = 0; i < DB_MAX; i++)
    {
        double SNR_db  = i;
        double SNR_num = pow(10, SNR_db / 10);
        double Eb      = 1;
        double N0      = Eb / SNR_num;
        double var     = N0 / 2;
        sigma[i]       = sqrt(var);
    }

    /*加噪声和衰落*/
    for (int i = 0; i < DB_MAX; i++) {
        ulong error = 0;
        for (int j = 0; j < N; j++) {
            int send        = gen_binomial_random(0.5);
            int in          = (send == 1) ? 1 : -1;
            double rayleigh = (fading == true) ? gen_rayleigh_random(1) : 1;
            double noise    = gen_normal_random(0, sigma[i]);
            double out      = rayleigh * in + noise;
            double recv     = (out > 0) ? 1 : 0;
            if (recv != send)
                error++;
        }
        snr_ber_p[i].SNR_db = i;
        snr_ber_p[i].BER    = (double)error / N;
    }
}

/*
 *产生0,1比特
 */
Serial* gen_signal(ulong N)
{
    srand(time(NULL));
    Serial* serial = (Serial*)malloc(sizeof(Serial));
    serial->N = N;
    serial->signal = (double*)malloc(sizeof(double) * N);
    for (ulong i = 0; i < N; i++)
        serial->signal[i] = gen_binomial_random(0.5);
    return serial;
}

void free_serial(Serial * serial)
{
    free(serial->signal);
    free(serial);
}

void add_noise(Serial * serial, double mean, double sigma)
{
    for (ulong i = 0; i < serial->N; i++)
        serial->signal[i] += gen_normal_random(mean, sigma);
}

void add_noise_and_fading(Serial * serial, double n_mean,
        double n_sigma, double r_sigma)
{
    for (ulong i = 0; i < serial->N; i++)
        serial->signal[i] = serial->signal[i] * gen_rayleigh_random(r_sigma)
            + gen_normal_random(n_mean, n_sigma);
}

/*串并变换, 变换为M路*/
Parallel* serial_to_parallel(Serial * serial, ushort M)
{
    Parallel * parallel = (Parallel*)malloc(sizeof(Parallel));
    //补零成M的倍数
    ulong old_N = serial->N;

    serial->N = serial->N + M - serial->N % M;
    serial->signal = (double*)realloc(serial->signal, sizeof(double) * serial->N);
    for (ulong n = old_N; n < serial->N; n++)
        serial->signal[n] = 0;

    parallel->M = M; parallel->N = serial->N / parallel->M;
    parallel->signal = (double**)malloc(sizeof(void*) * parallel->N);
    for (ulong n = 0; n < parallel->N; n++)
        parallel->signal[n] = (double*)malloc(sizeof(double) * parallel->M);

    for (ulong n = 0; n < parallel->N; n++) {
        for (ushort m = 0; m < parallel->M; m++)
            parallel->signal[n][m] = serial->signal[n * parallel->M + m];
    }

    return parallel;
}

void free_parallel(Parallel * parallel)
{
    for (ulong n = 0; n < parallel->N; n++)
        free(parallel->signal[n]);
    free(parallel->signal);
    free(parallel);
}

Serial* parallel_to_serial(Parallel * parallel)
{
    Serial* serial = (Serial*)malloc(sizeof(Serial));
    serial->N = parallel->N;
    serial->signal = (double*)malloc(sizeof(double) * serial->N);
    for (ulong n = 0; n < parallel->N; n++) {
        serial->signal[n] = 0;
        for (ushort m = 0; m < parallel->M; m++)
            serial->signal[n] += parallel->signal[n][m] * pow(2, m);
    }
    return serial;
}

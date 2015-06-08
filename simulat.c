#include "random.h"
#include "simulat.h"

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/*先计算噪声标准差*/
void init_sigma_array(double** sigma, int DB_MAX)
{
    *sigma = (double*)malloc(DB_MAX * sizeof(double));
    for (int i = 0; i < DB_MAX; i++)
    {
        double SNR_db  = i;
        double SNR_num = pow(10, SNR_db / 10);
        double Eb      = 1;
        double N0      = Eb / SNR_num;
        double var     = N0 / 2;
        (*sigma)[i]       = sqrt(var);
    }
}

void free_sigma_array(double* sigma)
{
    free(sigma);
}

void simulat_BPSK(const unsigned long N, const unsigned DB_MAX,
        const bool fading, const double* sigma, Point* snr_ber_p)
{
    /*加噪声和衰落*/
    for (int i = 0; i < DB_MAX; i++) {
        unsigned long error = 0;
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

void simulat_QPSK(const unsigned long N, const unsigned DB_MAX,
        const bool fading, const double* sigma, Point* snr_ber_p)
{
    for (int i = 0; i < DB_MAX; i++) {
        unsigned long error = 0;
        for (int I = 0, Q = 1; I < N && Q < N; I = I + 2, Q =Q + 2) {
            int send_I = gen_binomial_random(0.5);
            int send_Q = gen_binomial_random(0.5);
            int in_I = (send_I == 1) ? 1 : -1;
            int in_Q = (send_Q == 1) ? 1 : -1;
            double rayleigh_I = (fading == true) ? gen_rayleigh_random(1) : 1;
            double rayleigh_Q = (fading == true) ? gen_rayleigh_random(1) : 1;
            double noise_I = gen_normal_random(0, sigma[i]);
            double noise_Q = gen_normal_random(0, sigma[i]);
            double out_I = rayleigh_I * in_I + noise_I;
            double out_Q = rayleigh_Q * in_Q + noise_Q;
            double recv_I = (out_I > 0) ? 1 : 0;
            double recv_Q = (out_Q > 0) ? 1 : 0;
            if (recv_I != send_I || recv_Q != send_Q)
                error++;
        }
        snr_ber_p[i].SNR_db = i;
        snr_ber_p[i].BER    = (double)error / N;
    }
}

#include "random.h"
#include "simulat.h"

#include <stdbool.h>
#include <math.h>

void simulat_BPSK(const unsigned long N, const unsigned DB_MAX,
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

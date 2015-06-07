#include "random.h"
#include "simulat.h"

#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

void simulat_BPSK(const unsigned long N, const unsigned DB_MAX,
        const bool fading, SNR_BER* SNR_BER_p)
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
        SNR_BER_p[i].SNR_db = i;
        SNR_BER_p[i].BER    = (double)error / N;
    }
}

int main(void)
{
    srand(time(NULL));

    const unsigned long N = 10000000;
    const unsigned DB_MAX = 41;
    const bool fading     = true;
    double  sigma[DB_MAX];
    SNR_BER snr_ber[DB_MAX];
    simulat_BPSK(N, DB_MAX, fading, snr_ber);
    return 0;
}

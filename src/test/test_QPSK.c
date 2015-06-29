#include "simulat.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void simulat_QPSK(FILE* fp, ulong N, ulong M, int db, bool fading)
{
    srand(time(NULL));
    int m = 16;
    double sigma = sqrt( 1 / ( 2 * pow( 10, (double)db/10) ) );
    printf("==========伪随机序列==========\n");
    Serial* serial = gen_signal(N);
    for (ulong n = 0; n < serial->N; n++) {
        printf("%g\t", serial->signal[n]);
        if ((n + 1) % m == 0)
            printf("\n");
    }
    printf("\n");

    printf("==========QPSK==========\n");
    Serial_C* serial_c = QPSK(serial);
    for (ulong n= 0; n < serial_c->N; n++)
        printf("%g + %gj\n", serial_c->signal[n][0], serial_c->signal[n][1]);

    printf("==========加高斯噪声(和衰落)==========\n");
    if (fading)
        add_noise_and_fading_c(serial_c, 0, sigma, 1);
    else
        add_noise_c(serial_c, 0, sigma);
    for (ulong n= 0; n < serial_c->N; n++)
        printf("%g + %gj\n", serial_c->signal[n][0], serial_c->signal[n][1]);

    printf("==========抽样判决==========\n");
    judge4QPSK(serial_c);
    for (ulong n= 0; n < serial_c->N; n++)
        printf("%g + %gj\n", serial_c->signal[n][0], serial_c->signal[n][1]);

    printf("==========rQPSK==========\n");
    Serial* serial_1 = rQPSK(serial_c);
    for (ulong n = 0; n < serial_1->N; n++) {
        printf("%g\t", serial_1->signal[n]);
        if ((n + 1) % m == 0)
            printf("\n");
    }
    printf("\n");


    printf("==========统计误码率==========\n");
    ulong e = 0;
    for (ulong n = 0; n < serial->N; n++)
        if (serial->signal[n] != serial_1->signal[n])
            e++;
    fprintf(fp, "%d,%0.7f\n", db, (double)e/N);

    free(serial_1->signal);
    free(serial_1);
    free(serial_c->signal);
    free(serial_c);
    free(serial->signal);
    free(serial);
}

int main(void)
{
    FILE* fp;
    const ulong N = 1000000;
    const ulong M = 16;
    const int DB_MAX = 41;

    if ((fp = fopen("qpsk_no_fading.csv", "w+")) != NULL)  {
        fprintf(fp, "SNR_db, BER\n");
        for (int db = 0; db < DB_MAX; db++)
            simulat_QPSK(fp, N, M, db, false);
    }
    fclose(fp);

    if ((fp = fopen("qpsk_fading.csv", "w+")) != NULL) {
        fprintf(fp, "SNR_db, BER\n");
        for (int db = 0; db < DB_MAX; db++)
            simulat_QPSK(fp, N, M, db, true);
    }
    fclose(fp);

    return 0;
}

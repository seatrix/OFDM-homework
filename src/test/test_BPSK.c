#include "simulat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void simulat_BPSK(FILE* fp, ulong N, ulong M, int db, bool fading)
{
    srand(time(NULL));
    int m = 16;
    double sigma = sqrt( 1 / ( 2 * pow(10, (double)db/10) ) );
    printf("==========伪随机序列==========\n");
    Serial* serial = gen_signal(N);
    for (ulong n = 0; n < serial->N; n++) {
        printf("%g\t", serial->signal[n]);
        if ((n + 1) % m == 0)
            printf("\n");
    }
    printf("\n");

    printf("==========BPSK映射==========\n");
    Serial* serial_1 = BPSK(serial);
    for (ulong n = 0; n < serial_1->N; n++) {
        printf("%g\t", serial_1->signal[n]);
        if ((n + 1) % m == 0)
            printf("\n");
    }
    printf("\n");


    printf("==========加高斯噪声(和衰落)==========\n");
    if (fading)
        add_noise_and_fading(serial_1, 0, sigma, 0.5);
    else
        add_noise(serial_1, 0, sigma);

    for (ulong n = 0; n < serial_1->N; n++) {
        printf("%g\t", serial_1->signal[n]);
        if ((n + 1) % m == 0)
            printf("\n");
    }
    printf("\n");


    printf("==========抽样判决==========\n");
    judge4BPSK(serial_1);
    for (ulong n = 0; n < serial_1->N; n++) {
        printf("%g\t", serial_1->signal[n]);
        if ((n + 1) % m == 0)
            printf("\n");
    }
    printf("\n");


    printf("==========BPSK逆映射==========\n");
    rBPSK(serial_1);
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
    fprintf(fp, "%d,%lf\n", db, (double)e/N);


    free(serial_1->signal);
    free(serial_1);
    free(serial->signal);
    free(serial);
}

int main(void)
{
    FILE* fp;
    const ulong N = 10000000;
    const ulong M = 16;
    const int DB_MAX = 41;

    if ((fp = fopen("bpsk_no_fading.csv", "w+")) != NULL)  {
        fprintf(fp, "SNR_db, BER\n");
        for (int db = 0; db < DB_MAX; db++)
            simulat_BPSK(fp, N, M, db, false);
    }
    fclose(fp);

    if ((fp = fopen("bpsk_fading.csv", "w+")) != NULL)  {
        fprintf(fp, "SNR_db, BER\n");
        for (int db = 0; db < DB_MAX; db++)
            simulat_BPSK(fp, N, M, db, true);
    }
    fclose(fp);

    return 0;
}

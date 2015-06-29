#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "simulat.h"

void simulat_OFDM(FILE* fp, ulong N, ulong M, int db, bool fading)
{
    srand(time(NULL));
    double sigma = sqrt( 1 / ( 2 * ( pow( 10, (double)db/10 ) ) ) );

    printf("==========伪随机序列==========\n");
    Serial* serial = gen_signal(N);
    for (ulong n = 0; n < serial->N; n++) {
        printf("%g\t", serial->signal[n]);
        if ((n + 1) % M == 0)
            printf("\n");
    }
    printf("\n");
    printf("==========16QAM==========\n");
    Serial_C * serial_c = QAM16(serial);
    for (ulong n = 0; n < serial_c->N; n++)
        printf("%g + %gj\n", serial_c->signal[n][0], serial_c->signal[n][1]);

    printf("==========串并变换==========\n");
    Parallel_C* parallel_c = serial_to_parallel_c(serial_c, M);
    for (ulong n = 0; n < parallel_c->N; n++) {
        for (ulong m = 0; m < parallel_c->M; m++)
            printf("%g + %gj,", parallel_c->signal[n][m][0], parallel_c->signal[n][m][1]);
        printf("\n\n");
    }

    printf("==========ifft==========\n");
    ifft(parallel_c);
    for (ulong n = 0; n < parallel_c->N; n++) {
        for (ulong m = 0; m < parallel_c->M; m++)
            printf("%g + %gj\n", parallel_c->signal[n][m][0], parallel_c->signal[n][m][1]);
        printf("\n\n");
    }

    printf("==========加循环前缀==========\n");
    add_cycle_prefix(parallel_c);
    for (ulong n = 0; n < parallel_c->N; n++) {
        for (ulong m = 0; m < parallel_c->M; m++)
            printf("%g + %gj\n", parallel_c->signal[n][m][0], parallel_c->signal[n][m][1]);
        printf("\n\n");
    }

    printf("==========并串变换==========\n");
    Serial_C* serial_c_1 = parallel_to_serial_c(parallel_c);
    for (ulong n= 0; n < serial_c_1->N; n++)
        printf("%g + %gj\n", serial_c_1->signal[n][0], serial_c_1->signal[n][1]);

    printf("==========加高斯噪声(和衰落)==========\n");
    if (fading)
        add_noise_and_fading_c(serial_c, 0, sigma, 1);
    else
        add_noise_c(serial_c_1, 0, sigma);
    for (ulong n= 0; n < serial_c_1->N; n++)
        printf("%g + %gj\n", serial_c_1->signal[n][0], serial_c_1->signal[n][1]);

    printf("==========串并变换==========\n");
    Parallel_C* parallel_c_1 = serial_to_parallel_c(serial_c_1, M);
    for (ulong n = 0; n < parallel_c_1->N; n++) {
        for (ulong m = 0; m < parallel_c_1->M; m++)
            printf("%g + %gj,", parallel_c_1->signal[n][m][0], parallel_c_1->signal[n][m][1]);
        printf("\n\n");
    }

    printf("==========去循环前缀==========\n");
    del_cycle_prefix(parallel_c_1);
    for (ulong n = 0; n < parallel_c->N; n++) {
        for (ulong m = 0; m < parallel_c->M; m++)
            printf("%g + %gj\n", parallel_c->signal[n][m][0], parallel_c->signal[n][m][1]);
        printf("\n\n");
    }

    printf("==========fft==========\n");
    fft(parallel_c_1);
    for (ulong n = 0; n < parallel_c->N; n++) {
        for (ulong m = 0; m < parallel_c->M; m++)
            printf("%g + %gj\n", parallel_c->signal[n][m][0], parallel_c->signal[n][m][1]);
        printf("\n\n");
    }

    printf("==========并串变换==========\n");
    Serial_C* serial_c_2 = parallel_to_serial_c(parallel_c);
    for (ulong n = 0; n < serial_c_2->N; n++)
        printf("%g + %gj\n", serial_c_2->signal[n][0], serial_c_2->signal[n][1]);

    printf("==========抽样判决==========\n");
    judge4QAM16(serial_c_2);
    for (ulong n = 0; n < serial_c_2->N; n++)
        printf("%g + %gj\n", serial_c_2->signal[n][0], serial_c_2->signal[n][1]);

    printf("==========r16QAM==========\n");
    Serial* serial_3 = rQAM16(serial_c);
    for (ulong n = 0; n < N; n++) {
        printf("%g\t", serial_3->signal[n]);
        if ((n + 1) % M == 0)
            printf("\n");
    }
    printf("\n");

    printf("==========统计误码率==========\n");
    ulong e = 0;
    for (ulong n = 0; n < serial->N; n++)
        if (serial->signal[n] != serial_3->signal[n])
            e++;
    fprintf(fp, "%d,%0.7f\n", db, (double)e/N);


    free(serial_3->signal);
    free(serial_3);
    free(serial_c_2);
    free(parallel_c_1->signal);
    free(parallel_c_1);
    free(serial_c_1);
    free(parallel_c->signal);
    free(parallel_c);
    free(serial_c->signal);
    free(serial_c);
    free(serial->signal);
    free(serial);
}

int main(void)
{
    FILE *fp;
    const ulong N = 1000000;
    const ulong M = 16;
    const int DB_MAX = 41;
    if ((fp = fopen("ofdm_no_fading.csv", "w+")) != NULL) {
        fprintf(fp, "SNR_db, BER\n");
        for (int db = 0; db < DB_MAX; db++)
            simulat_OFDM(fp, N, M, db, false);
    }
    fclose(fp);

    if ((fp = fopen("ofdm_fading.csv", "w+")) != NULL)  {
        fprintf(fp, "SNR_db, BER\n");
        for (int db = 0; db < DB_MAX; db++)
            simulat_OFDM(fp, N, M, db, true);
    }
    fclose(fp);

    return 0;
}

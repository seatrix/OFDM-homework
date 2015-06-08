#include "simulat.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

/*
 *开始仿真, 结果保存文件
 */
void simulat(int DB_MAX, unsigned long N)
{
    FILE *fp_BPSK_no_fading = fopen("./BPSK_no_fading.csv", "w+");//加噪声的结果文件
    FILE *fp_BPSK_fading    = fopen("./BPSK_fading.csv", "w+");//加衰落和噪声的结果文件
    FILE *fp_QPSK_no_fading = fopen("./QPSK_no_fading.csv", "w+");
    FILE *fp_QPSK_fading    = fopen("./QPSK_fading.csv", "w+");
    if (fp_BPSK_no_fading == NULL ||
        fp_BPSK_fading    == NULL ||
        fp_QPSK_no_fading == NULL ||
        fp_QPSK_fading    == NULL) {
        printf("文件创建失败\n");
        return;
    }
    srand(time(NULL));
    Point snr_ber[DB_MAX];
    double* sigma = NULL;
    init_sigma_array(&sigma, DB_MAX);
    /*初始化工作结束*/

    simulat_BPSK(N, DB_MAX, true, sigma, snr_ber);
    fprintf(fp_BPSK_fading, "SNR_db,BER\n");
    for (int i = 0; i < DB_MAX; i++)
        fprintf(fp_BPSK_fading, "%f,%f\n", snr_ber[i].SNR_db, snr_ber[i].BER);

    simulat_BPSK(N, DB_MAX, false, sigma, snr_ber);
    fprintf(fp_BPSK_no_fading, "SNR_db,BER\n");
    for (int i = 0; i < DB_MAX; i++)
        fprintf(fp_BPSK_no_fading, "%f,%f\n", snr_ber[i].SNR_db, snr_ber[i].BER);

    simulat_QPSK(N, DB_MAX, true, sigma, snr_ber);
    fprintf(fp_QPSK_fading, "SNR_db,BER\n");
    for (int i = 0; i < DB_MAX; i++)
        fprintf(fp_QPSK_fading, "%f,%f\n", snr_ber[i].SNR_db, snr_ber[i].BER);

    simulat_QPSK(N, DB_MAX, false, sigma, snr_ber);
    fprintf(fp_QPSK_no_fading, "SNR_db,BER\n");
    for (int i = 0; i < DB_MAX; i++)
        fprintf(fp_QPSK_no_fading, "%f,%f\n", snr_ber[i].SNR_db, snr_ber[i].BER);

    /*收尾工作开始*/
    free_sigma_array(sigma);
    fclose(fp_BPSK_fading);
    fclose(fp_BPSK_no_fading);
    fclose(fp_QPSK_fading);
    fclose(fp_QPSK_no_fading);
}

int main(void)
{
    const int DB_MAX = 41;
    const unsigned long N = 100000;
    simulat(DB_MAX, N);
    return 0;
}


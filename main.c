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
    FILE *fp_no_fading = fopen("./no_fading.csv", "w+");//加噪声的结果文件
    FILE *fp_fading = fopen("./fading.csv", "w+");//加衰落和噪声的结果文件
    if (fp_no_fading == NULL)
    {
        printf("加噪声文件创建失败\n");
        return;
    }
    if (fp_fading == NULL)
    {
        printf("加衰落和噪声文件创建失败\n");
    }

    srand(time(NULL));

    Point snr_ber[DB_MAX];

    simulat_BPSK(N, DB_MAX, true, snr_ber);
    fprintf(fp_fading, "SNR_db,BER\n");
    for (int i = 0; i < DB_MAX; i++)
        fprintf(fp_fading, "%f,%f\n", snr_ber[i].SNR_db, snr_ber[i].BER);

    simulat_BPSK(N, DB_MAX, false, snr_ber);
    fprintf(fp_no_fading, "SNR_db,BER\n");
    for (int i = 0; i < DB_MAX; i++)
        fprintf(fp_no_fading, "%f,%f\n", snr_ber[i].SNR_db, snr_ber[i].BER);

    fclose(fp_fading);
    fclose(fp_no_fading);
}

int main(void)
{
    const int DB_MAX = 41;
    const unsigned long N = 10000000;
    simulat(DB_MAX, N);
    return 0;
}


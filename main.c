#include "librandom.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>


void run(FILE *fp, int DB_MAX, unsigned long N, bool fading)
{
    fprintf(fp, "SNR_db,BER\n");//csv文件的头部
    srand(time(NULL));

    /*
     *信噪比从0db到DB_MAX,
     *信号功率恒定且为1w,
     *噪声功率为1/SRN_num, 即噪声方差.
     */
    double SNR_db[DB_MAX];
    double SNR_num[DB_MAX];
    double var[DB_MAX];
    for (int i = 0; i < DB_MAX; i++)
    {
        SNR_db[i] = i;
        SNR_num[i] = pow(10, SNR_db[i] / 10);
        var[i] = 1 / SNR_num[i];
    }

    double BER[DB_MAX];
    for (int i = 0; i < DB_MAX; i++)
    {
        unsigned long error = 0;
        for (int j = 0; j < N; j++)
        {
            int b = gen_binomial_random(0.5);//0,1比特
            int in = (b == 1) ? 1 : -1;//输入
            double r = (fading == true) ? gen_rayleigh_random(1) : 1;//瑞利衰落
            double n = gen_normal_random(0, sqrt(var[i]));//高斯噪声
            double out = r * in + n;//输出
            double v = (out > 0) ? 1 : 0;//判决
            if (v != b)
                error++;
        }
        BER[i] = (double)error / N;
        fprintf(fp, "%d,%lf\n", i, BER[i]);
    }
}

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

    run(fp_fading, DB_MAX, N, true);
    run(fp_no_fading, DB_MAX, N, false);


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


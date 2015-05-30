#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(void)
{
    srand(time(NULL));

    /*
     *信噪比从0db到40db,
     *信号功率恒定且为1w,
     *噪声功率为1/Eb2N_num, 即噪声方差.
     */
    const int DB_MAX = 40;
    double Eb2N_db[DB_MAX];
    double Eb2N_num[DB_MAX];
    double var[DB_MAX];
    for (int i = 0; i < DB_MAX; i++)
    {
        Eb2N_db[i] = i;
        Eb2N_num[i] = pow(10, Eb2N_db[i] / 10);
        var[i] = 1 / Eb2N_num[i];
    }

    const int N = 10000000;
    double BER[DB_MAX];
    for (int i = 0; i < DB_MAX; i++)
    {
        int error = 0;
        for (int j = 0; j < N; j++)
        {
            int b = gen_binomial_random(0.5);//0,1比特
            int in = (b == 1) ? 1 : -1;//输入
            double r = gen_rayleigh_random();//瑞利衰落
            double n = gen_normal_random(0, sqrt(var[i]));//高斯噪声
            double out = r * in + n;//输出
            double v = (out > 0) ? 1: 0;//判决
            if (v != b)
                error++;
        }
        BER[i] = (double)error / N;
        printf("S2N: %d, BER: %lf\n", i, BER[i]);
    }

    return 0;
}


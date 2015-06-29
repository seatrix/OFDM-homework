#include "simulat.h"
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

/*16QAM星座图*/
static const double _0[4]  = {0, 0, 0, 0};
static const double _1[4]  = {0, 0, 0, 1};
static const double _2[4]  = {0, 0, 1, 0};
static const double _3[4]  = {0, 0, 1, 1};
static const double _4[4]  = {0, 1, 0, 0};
static const double _5[4]  = {0, 1, 0, 1};
static const double _6[4]  = {0, 1, 1, 0};
static const double _7[4]  = {0, 1, 1, 1};
static const double _8[4]  = {1, 0, 0, 0};
static const double _9[4]  = {1, 0, 0, 1};
static const double _10[4] = {1, 0, 1, 0};
static const double _11[4] = {1, 0, 1, 1};
static const double _12[4] = {1, 1, 0, 0};
static const double _13[4] = {1, 1, 0, 1};
static const double _14[4] = {1, 1, 1, 0};
static const double _15[4] = {1, 1, 1, 1};

/*每个bit的能量是1, 则每个符号的能量是4, 采用方形16QAM, 则振幅分别可以选为1/3 2/3 1*/
#define _x1 1
#define _x2 1.7320508075688772
#define _y1 1
#define _y2 1.7320508075688772

static const double __0[2]  = {_x1,   _y1};
static const double __1[2]  = {_x2,   _y1};
static const double __2[2]  = {_x1,   _y2};
static const double __3[2]  = {_x2,   _y2};
static const double __4[2]  = {_x1,  -_y1};
static const double __5[2]  = {_x1,  -_y2};
static const double __6[2]  = {_x2,  -_y1};
static const double __7[2]  = {_x2,  -_y2};
static const double __8[2]  = {-_x1,  _y1};
static const double __9[2]  = {-_x1,  _y2};
static const double __10[2] = {-_x2,  _y1};
static const double __11[2] = {-_x2,  _y2};
static const double __12[2] = {-_x1, -_y1};
static const double __13[2] = {-_x2, -_y1};
static const double __14[2] = {-_x1, -_y2};
static const double __15[2] = {-_x2, -_y2};


void simulat_BPSK(const ulong N, const ulong DB_MAX,
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
        ulong error = 0;
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

/*
 *产生0,1比特
 */
Serial* gen_signal(ulong N)
{
    srand(time(NULL));
    Serial* serial = (Serial*)malloc(sizeof(Serial));
    serial->N = N;
    serial->signal = (double*)malloc(sizeof(double) * N);
    for (ulong i = 0; i < N; i++)
        serial->signal[i] = gen_binomial_random(0.5);
    return serial;
}


void bipolar(Serial * serial, double scale)
{
    for (ulong n = 0; n < serial->N; n++) {
        if (serial->signal[n] == 1)
            serial->signal[n] = scale;
        else
            serial->signal[n] = -scale;
    }


}

void add_noise(Serial * serial, double mean, double sigma)
{
    for (ulong n = 0; n < serial->N; n++)
        serial->signal[n] += gen_normal_random(mean, sigma);
}

void add_noise_c(Serial_C * serial_c, double mean, double sigma)
{
    for (ulong n = 0; n < serial_c->N; n++) {
        serial_c->signal[n][0] += gen_normal_random(mean, sigma);
        serial_c->signal[n][1] += gen_normal_random(mean, sigma);
    }
}

void add_noise_and_fading(Serial * serial, double n_mean,
        double n_sigma, double r_sigma)
{
    for (ulong n = 0; n < serial->N; n++)
        serial->signal[n] = serial->signal[n] * gen_rayleigh_random(r_sigma)
            + gen_normal_random(n_mean, n_sigma);
}

void add_noise_and_fading_c(Serial_C * serial_c, double n_mean,
        double n_sigma, double r_sigma)
{
    for (ulong n = 0; n < serial_c->N; n++) {
        serial_c->signal[n][0] = serial_c->signal[n][0] * gen_rayleigh_random(r_sigma)
            + gen_normal_random(n_mean, n_sigma);
        serial_c->signal[n][1] = serial_c->signal[n][1] * gen_rayleigh_random(r_sigma)
            + gen_normal_random(n_mean, n_sigma);
    }
}

/*串并变换, 变换为M路*/
Parallel* serial_to_parallel(Serial * serial, ulong M)
{
    Parallel * parallel = (Parallel*)malloc(sizeof(Parallel));
    //补零成M的倍数
    ulong old_N = serial->N;
    serial->N = (serial->N % M == 0) ? serial->N : serial->N + M - serial->N % M;
    serial->signal = (double*)realloc(serial->signal, sizeof(double) * serial->N);
    for (ulong n = old_N; n < serial->N; n++)
        serial->signal[n] = 0;

    parallel->M = M; parallel->N = serial->N / parallel->M;
    parallel->signal = (double**)malloc(sizeof(double*) * parallel->N);
    for (ulong n = 0; n < parallel->N; n++)
        parallel->signal[n] = &serial->signal[n * parallel->M];
    return parallel;
}

/*串并变换, 变换为M路, 复数版本*/
Parallel_C* serial_to_parallel_c(Serial_C* serial_c, ulong M)
{
    Parallel_C* parallel_c = (Parallel_C*)malloc(sizeof(Parallel_C));
    //补零成M的倍数
    ulong old_N = serial_c->N;
    serial_c->N = (serial_c->N % M == 0) ? serial_c->N : serial_c->N + M - serial_c->N % M;
    serial_c->signal = (fftw_complex*)realloc(serial_c->signal, sizeof(fftw_complex) * serial_c->N);
    for (ulong n = old_N; n < serial_c->N; n++) {
        serial_c->signal[n][0] = 0;
        serial_c->signal[n][1] = 0;
    }
    parallel_c->M = M; parallel_c->N = serial_c->N / parallel_c->M;
    parallel_c->signal = (fftw_complex**)malloc(sizeof(fftw_complex*) * parallel_c->N);
    for (ulong n = 0; n < parallel_c->N; n++)
        parallel_c->signal[n] = &serial_c->signal[n * parallel_c->M];
    return parallel_c;
}

void add_cycle_prefix(Parallel_C* parallel_c)
{
    /*
    设计时使用的一维连续数组代表二维数组, 现在无法改变子数组宽度, 先注释掉
    ulong nsize = parallel_c->M * (1 + 0.25);
    for (ulong n = 0; n < parallel_c->N; n++) {
        parallel_c->signal[n] = (fftw_complex*)realloc(parallel_c->signal[n], sizeof(fftw_complex) * nsize);
        for (ulong m = parallel_c->M; m < nsize; m++) {
            parallel_c->signal[n][m][0] = parallel_c->signal[n][m - nsize][0];
            parallel_c->signal[n][m][1] = parallel_c->signal[n][m - nsize][0];
        }
    }
    parallel_c->M = nsize;
    */
}

void del_cycle_prefix(Parallel_C* parallel_c)
{
    /*
    ulong nsize = parallel_c->M / (1 + 0.25);
    for (ulong n = 0; n < parallel_c->N; n++)
        parallel_c->signal[n] = (fftw_complex*)realloc(parallel_c->signal[n], sizeof(fftw_complex) * nsize);
    parallel_c->M = nsize;
    */
}

Serial_C* parallel_to_serial_c(Parallel_C* parallel_c)
{
    Serial_C* serial_c = (Serial_C*)malloc(sizeof(Serial_C));
    serial_c->N = parallel_c->N * parallel_c->M;
    serial_c->signal = *parallel_c->signal;
    return serial_c;
}


/*每个bit的能量是1, 则每个符号的能量是2, 所以模是sqrt(2)*/
Serial_C* QPSK(Serial* serial)
{
    Serial_C* serial_c = (Serial_C*)malloc(sizeof(Serial_C));
    serial_c->N = (serial_c->N % 2 == 0) ? serial_c->N : serial->N + 2 - serial->N % 2;
    serial_c->signal = (fftw_complex*)malloc(sizeof(fftw_complex) * serial_c->N);
    for (ulong n = 0; n < serial_c->N; n++)
        for (ulong m = 0; m < 2; m++) {
            if (n * 2 + m < serial->N)
                serial_c->signal[n][m] = serial->signal[n * 2 + m];
            else
                serial_c->signal[n][m] = 0;
        }

    return serial_c;
}

Serial* rQPSK(Serial_C* serial_c)
{
    Serial* serial = (Serial*)malloc(sizeof(Serial));
    serial->N = serial_c->N * 2;
    serial->signal = (double*)malloc(sizeof(double) * serial->N);
    for (ulong n = 0; n < serial_c->N; n++) {
        serial->signal[2 * n] = serial_c->signal[n][0];
        serial->signal[2 * n + 1] = serial_c->signal[n][1];
    }
    return serial;
}

Serial_C* QAM16(Serial* serial)
{
    Parallel* parallel4 = serial_to_parallel(serial, 4);

    Serial_C* serial_c = (Serial_C*)malloc(sizeof(Serial_C));
    serial_c->N = parallel4->N;
    serial_c->signal = (fftw_complex*)malloc(sizeof(fftw_complex) * serial_c->N);

    for (ulong n = 0; n < serial_c->N; n++) {
        /*fuck C, 不能switch 数组!*/
        if (memcmp(parallel4->signal[n], _0 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __0, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _1 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __1, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _2 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __2, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _3 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __3, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _4 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __4, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _5 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __5, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _6 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __6, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _7 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __7, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _8 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __8, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _9 , sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __9, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _10, sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __10, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _11, sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __11, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _12, sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __12, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _13, sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __13, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _14, sizeof(double) * 4) == 0) 
            memcpy(serial_c->signal[n], __14, sizeof(double) * 2);
        else if (memcmp(parallel4->signal[n], _15, sizeof(double) * 4) == 0)
            memcpy(serial_c->signal[n], __15, sizeof(double) * 2);
        else 
            exit(-1);
    }
    free(parallel4->signal);
    free(parallel4);
    return serial_c;
}

/*离谁的距离最小*/
static double min_d(double a, double es[], int N)
{
    double *d = (double*)malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) d[i] = fabs(a - es[i]);
    double dmin = d[0]; int j = 0;
    for (int i = 1; i < N; i++)
        if (d[i] < dmin) { dmin = d[i]; j = i; }
    free(d);
    return j;
}

/*判决*/
void judge4QAM16(Serial_C* serial_c)
{
    double es[4] = {-_x2, -_x1, _x1, _x2};
    for (ulong n = 0; n < serial_c->N; n++) {
        serial_c->signal[n][0] = min_d(serial_c->signal[n][0], es, 4);
        serial_c->signal[n][1] = min_d(serial_c->signal[n][1], es, 4);
    }
}

Serial* rQAM16(Serial_C* serial_c)
{
    judge4QAM16(serial_c);
    Serial* serial = (Serial*)malloc(sizeof(Serial));
    serial->N = serial_c->N * 2;
    serial->signal = (double*)malloc(sizeof(double) * serial->N);
    for (ulong n =  0; n < serial_c->N; n++) {
        if (memcmp(serial_c->signal[n], __0, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _0, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __1, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _1, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __2, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _2, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __3, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _3, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __4, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _4, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __5, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _5, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __6, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _6, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __7, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _7, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __8, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _8, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __9, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _9, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __10, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _10, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __11, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _11, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __12, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _12, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __13, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _13, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __14, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _14, sizeof(double) * 4);
        else if (memcmp(serial_c->signal[n], __15, sizeof(double) * 2))
            memcpy(&serial->signal[4 * n], _15, sizeof(double) * 4);
        else
            exit(-1);
    }
    return serial;
}

static void transform(Parallel_C* parallel_c, int forward)
{
    fftw_complex in[parallel_c->M];
    fftw_complex out[parallel_c->M];
    fftw_plan p = fftw_plan_dft_1d(parallel_c->M, in, out, forward, FFTW_MEASURE);
    for (ulong n = 0; n < parallel_c->N; n++) {
        memcpy(in, parallel_c->signal[n], sizeof(fftw_complex) * parallel_c->M);
        fftw_execute(p);
        memcpy(parallel_c->signal[n], out, sizeof(fftw_complex) * parallel_c->M);
    }
    fftw_destroy_plan(p);
}

void ifft(Parallel_C* parallel_c)
{
    transform(parallel_c, 1);
    for (ulong n = 0; n < parallel_c->N; n++)
        for (ulong m = 0; m < parallel_c->M; m++) {
            parallel_c->signal[n][m][0] /= 16;
            parallel_c->signal[n][m][1] /= 16;
        }
}

void fft(Parallel_C* parallel_c)
{
    transform(parallel_c, -1);
}

static void swap(ulong *x, ulong *y)
{
    ulong t = *x;
    *x = *y;
    *y = t;
}

void transpose_array(Parallel* parallel)
{
    double tmparr[parallel->M][parallel->N];
    for (ulong n = 0; n < parallel->N; n++)
        for (ulong m = 0; m < parallel->M; m++)
            tmparr[m][n] = parallel->signal[n][m];
    memcpy(*(parallel->signal), &tmparr, sizeof(double) * parallel->N * parallel->M);

    swap(&(parallel->N), &(parallel->M));

    double** tmp = (double**)malloc(sizeof(double*) * parallel->N);
    for (ulong n = 0; n < parallel->N; n++)
        tmp[n] = &((*(parallel->signal))[n * parallel->M]);
    free(parallel->signal);
    parallel->signal = tmp;
}

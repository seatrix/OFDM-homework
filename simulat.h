#ifndef SIMULAT_H
#define SIMULAT_H

#include <stdbool.h>

typedef struct Point {
    double SNR_db;
    double BER;
} Point;


/*
 *参数:
 *    N:         产生的信号比特数
 *    DB_MAX:    SNR_db最大值
 *    fading:    是否添加瑞利衰落
 *    SNR_BER_p: 信噪比/误码率--作为返回值, 指向Point结构组成的数组
 */

void simulat_BPSK(const unsigned long N, const unsigned DB_MAX,
        const bool fading, const double* sigma, Point* snr_ber_p);


void simulat_QPSK(const unsigned long N, const unsigned DB_MAX,
        const bool fading, const double* sigma, Point* snr_ber_p);

/*注意参数是double**, 如果传入double*, 传入的是指针的副本, 
 * 函数体内改变指针指向的地址, 并不会反映到函数体外*/
void init_sigma_array(double** sigma, int DB_MAX);

void free_sigma_array(double* sigma);

#endif

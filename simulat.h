#ifndef SIMULAT_H
#define SIMULAT_H

#include <stdbool.h>

typedef struct SNR_BER {
    int SNR_db;
    int BER;
} SNR_BER;


/*
 *参数:
 *    N:         产生的信号比特数
 *    DB_MAX:    SNR_db最大值
 *    fading:    是否添加瑞利衰落
 *    SNR_BER_p: 信噪比/误码率--作为返回值, 指向SNR_BER结构组成的数组
 */

void simulat_BPSK(const unsigned long N, const unsigned DB_MAX,
        const bool fading, SNR_BER* SNR_BER_p);

#endif

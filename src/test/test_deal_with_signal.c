#include <stdio.h>
#include "simulat.h"

int main(void)
{
    ulong N = 60;
    ushort M = 16;
    Serial* serial = gen_signal(N);
    Parallel* parallel = serial_to_parallel(serial, M);
    for (ulong n = 0; n < parallel->N; n++) {

        for (ushort m = 0; m < parallel->M; m++) {
            if (parallel->signal[n][m] == serial->signal[n * parallel->M + m])
                printf("%g", parallel->signal[n][m]);
            else {
                printf("error!\n");
                goto error;
            }
        }
        printf("\n");
    }

    Serial* serial_baseM = parallel_to_serial(parallel);
    for (int i = 0; i < serial_baseM->N; i++)
        printf("%g ", serial_baseM->signal[i]);
    /*todo: 将输出的数与Python计算的进行比较*/

    free_serial(serial_baseM);
error:
    free_serial(serial);
    free_parallel(parallel);
    return 0;
}

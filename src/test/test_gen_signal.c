#include <stdio.h>
#include "simulat.h"

int main(void)
{
    ulong N = 100;
    ushort M = 2;
    Serial* serial = gen_signal(N);
    Parallel* parallel = serial_to_parallel(serial, M);
    for (ulong n = 0; n < parallel->N; n++)
        for (ushort m = 0; m < parallel->M; m++)
            if (parallel->signal[n][m] != serial->signal[n * parallel->M + m]) {
                printf("error!\n");
                goto error;
            }
error:
    free_serial(serial);
    free_parallel(parallel);
    return 0;
}

#include <stdlib.h>
#include <math.h>

/*
 *产生(0-1)均匀分布的随机数
 *http://stackoverflow.com/questions/11641629/generating-a-uniform-distribution-of-integers-in-c
 */
double gen_uniform_random(void)
{
    return rand() / (RAND_MAX + 1.);
}


/*
 *z = r * cos(phi)
 *w = r * sin(phi)
 *phi = pi * (2u - 1)
 *r = sqrt(-2ln(v))
 *u,v是(0-1)区间上均匀分布的随机数
 *r服从瑞利分布
 *z, w都服从标准正态分布
 */

/*
 *瑞利分布
 *f(x) = ( x / (a^2) ) * e^( -x^2 / (2*a^2) ) if x >=0 else 0
 *令 a=1
 */
double gen_rayleigh_random(double sigma)
{
    double v = gen_uniform_random();
    return sigma * sqrt(-2 * log(v));
}

/*
 *标准正态分布
 */
double gen_standard_normal_random(void)
{
    double u = gen_uniform_random();
    return gen_rayleigh_random(1) * cos(M_PI * ( 2 * u - 1 ));
}

/*
 *正态分布
 *如果X~N(0, 1)
 *那么aX+u ~ N(u, a^2)
 */
double gen_normal_random(double mean, double sigma)
{
    return sigma * gen_standard_normal_random() + mean;
}


/*
 *产生取值为0, 1的二项分布的随机数, 返回0的概率为p, 返回1的概率为1-p
 */
int gen_binomial_random(double p)
{
    return gen_uniform_random() > p ? 1: 0;
}


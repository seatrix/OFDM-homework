#ifndef LIBRANDOM_H
#define LIBRANDOM_H

double gen_uniform_random(void);
double gen_rayleigh_random(double sigma);
double gen_standard_normal_random(void);
double gen_normal_random(double mean, double sigma);
int gen_binomial_random(double p);

#endif

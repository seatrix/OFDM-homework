#ifndef RANDOM_H
#define RANDOM_H

double gen_uniform_random();
double gen_rayleigh_random();
double gen_standard_normal_random();
double gen_normal_random(double mean, double sigma);
int gen_binomial_random(double p);

#endif

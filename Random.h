/*
	Random.h
	Header file for methods using random numbers, specifically normally distributed numbers
*/

#ifndef RANDOM_H
#define RANDOM_H

#define PI 3.14159265358979323846264338327

#include "mt19937ar.h"

// Returns a random normally distributed number, mean 0, standard deviation 1
double box_muller_rand(){
	double r1 = genrand_real3();
	double r2 = genrand_real3();
	return sqrt(-2 * log(r1)) * cos(2 * PI * r2);
}

#endif
/*
	Random.h
	Header file for methods using random numbers, specifically normally distributed numbers
*/

#ifndef RANDOM_H
#define RANDOM_H

#define PI 3.14159265358979323846264338327


// Returns a uniformly distributed random number in range (0,1)
double uniform_rand(){
	return (((double)rand()) / RAND_MAX);
}

// Returns a random normally distributed number, mean 0, standard deviation 1
double box_muller_rand(){
	double r1 = uniform_rand();
	double r2 = uniform_rand();
	return sqrt(-2 * log(r1)) * cos(2 * PI * r2);
}

#endif
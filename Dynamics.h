/*
	Dynamics.h
	Header file for how the chosen dynamics changes the position (and possibly momentum)
	of a random walker.

	Reference for BAOAB method: 
	"Rational Construction of Stochastic Numerical Methods for Molecular Sampling" by
	Benedict Leimkuhler and Charles Matthews. arXiv: https://arxiv.org/abs/1203.5428
*/

#ifndef DYNAMICS_H
#define DYNAMICS_H
#include <math.h>
#include "Parameters.h"
#include "Random.h"
#include "MiscFunctions.h"
#include "mt19937ar.h"

#define PI 3.14159265358979323846264338327

// Typedef for a function pointer
typedef double (*DynamicsFun)(double, parameters*);

// High friction limit of BAOAB method
double BAOAB_limit(double x, parameters *params){
    // Set x according to the formula given in Leimkuhler & Matthews
	x = x - params->timestep * (*params->Poten_deriv)(x) / params->mass +\
		sqrt(0.5 * params->kT * params->timestep / params->mass) * (params->R[0] + params->R[1]);

	params->R[0] = params->R[1]; // Current value of R becomes the next one
	params->R[1] = box_muller_rand(); // Next value for R is drawn from a normal distribution
	return x;
}


// // Regular BAOAB method
// void BAOAB_regular(double ret_arr[2], double x, double p, parameters *params, double R[2]){
// 	double p_half = p - 0.5 * params->timestep * params->Poten_deriv(x);
// 	double x_half = x + 0.5 * params->timestep * p_half / params->mass;
// 	double p_half_bar = params->const1 * p_half + params->const3 * sqrt(params->mass) * R[1];
// 	ret_arr[0] = x_half + 0.5 * params->timestep * p_half_bar / params->mass;
// 	ret_arr[1] = p_half_bar - 0.5 * params->timestep * params->Poten_deriv(x);
// }


// Using Monte-Carlo method to determine where the particle jumps to
double Monte_Carlo_step(double x, parameters *params){
    // Possible new position, uniformly distributed
	double new_x = x + (2 * genrand_real1() - 1) * params->jump_size;

    // Difference in the potential between the new position and the current position
	double potential_difference = params->Poten(new_x) - params->Poten(x);

	double P_move = min(1, exp(-potential_difference / params->kT)); // Probability of moving to new_x
	if (genrand_real1() < P_move){
        // Moves to new_x with probability P_move
		return new_x;
	}
	return x; // If not moving to new_x, stay where currently are
}

DynamicsFun Dynamics_selector(char name[]){
    DynamicsFun dyn_fun; // Function for the chosen dynamics
	if (strcmp(name, "BAOAB_LIMIT") == 0){
		dyn_fun = &BAOAB_limit;
	}
	else if (strcmp(name, "MONTE-CARLO") == 0){
		dyn_fun = &Monte_Carlo_step;
	}
	return dyn_fun;
}


#endif

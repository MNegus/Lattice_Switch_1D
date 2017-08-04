#ifndef DYNAMICS_H
#define DYNAMICS_H
#include "Parameters.h"

double BAOAB_limit(double x, parameters *params, double R[2]){
	return x - params->timestep * (*params->Poten_deriv)(x) / params->mass +\
		sqrt(0.5 * params->kT * params->timestep / params->mass) * (R[0] + R[1]);
}

void BAOAB_regular(double ret_arr[2], double x, double p, parameters *params, double R[2]){
	double p_half = p - 0.5 * params->timestep * params->Poten_deriv(x);
	double x_half = x + 0.5 * params->timestep * p_half / params->mass;
	double p_half_bar = params->const1 * p_half + params->const3 * sqrt(params->mass) * R[1];
	ret_arr[0] = x_half + 0.5 * params->timestep * p_half_bar / params->mass;
	ret_arr[1] = p_half_bar - 0.5 * params->timestep * params->Poten_deriv(x);
}

#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "Potentials.h"

#define PI 3.14159265358979323846264338327

/* Structure to store relevant parameters */
struct parameters{
	// Numerical parameters
	double tot_timesteps;
	double timestep;
	double sample_portion;
	int start_well;
	int switch_regularity;
	
	// Physical parameters
	double kT;
	double mass;
	double shift_value;
	double minima[];
};
typedef struct parameters parameters;


// Returns the x-position of a particle in space given its displacement from the current well minima
double x_pos(parameters *params, int well, double displacement){
	return displacement + params->minima[well];
}

// Returns the displacement from the current well minima given its x position in space
double well_dis(parameters *params, int well, double x_position){
	return x_position - params->minima[well];
}

// Returns a uniformly distributed random number in range (0,1)
double uniform_rand(){
	return (((double)rand()) / RAND_MAX);
}

// Creates two normally distributed numbers, mean 0, standard deviation 1
void box_muller_rand(double ret_arr[]){
	double r1 = uniform_rand();
	double r2 = uniform_rand();
	ret_arr[0] = sqrt(-2 * log(r1)) * cos(2 * PI * r2);
	ret_arr[1] = sqrt(-2 * log(r1)) * sin(2 * PI * r2);
}

// Returns the minimum of two double values
double min(double x1, double x2){
	if (x1 < x2) {
		return x1;
	}
	else {
		return x2;
	}
}

double BAOAB_limit(double x, parameters *params, PotentialFun Poten_deriv, double R[]){
	return x - params->timestep * (*Poten_deriv)(x) / params->mass +\
		sqrt(0.5 * params->kT * params->timestep / params->mass) * (R[0] + R[1]);
}

// Calculates the free energy different between states in the two wells of a given potential function
void calc_energy_diff(double ret_arr[], char potential_name[], parameters *params, char data_dir, char output_filename){
	srand ( time(NULL) ); // Seeds the random number generators
	
	/* Retrieves constants and functions relevant to the specific potential chosen */
	double const_arr[] = {0, 0, 0}; // Array for constants
	PotentialFun func_arr[] = {0, 0, 0}; // Array for potential functions 
	U_selector(const_arr, func_arr, potential_name); // Fills the constant and function arrays

	params->minima[0] = const_arr[0]; // x-coordinates of the minima of the wells
	params->minima[1] = const_arr[1];
	params->shift_value = const_arr[2]; // The amount the right minima has been shifted upwards

	PotentialFun Poten = func_arr[0]; // Function pointer for the potential
	PotentialFun Poten_shift = func_arr[1]; // Function pointerfor the shifted potential
	PotentialFun Poten_deriv = func_arr[2]; // Function pointer for the derivative of the potential

	double x = x_pos(params, params->start_well, 0); // x-position initially at the bottom of the starting well

	double no_left = 0; // Number of timesteps that the particle is in the left well

	int cur_well = params->start_well; // Indicates which well the particle is in (0 is left well, 1 is right well)

	// Generates normally distributed values for R, which stores the current value and the value at the next timestep
	double normal_dist[] = {0, 0};
	box_muller_rand(normal_dist);
	double R[] = {normal_dist[0], normal_dist[1]};

	// Perform lattice switching method
	for (long m=0; m < params->tot_timesteps; m++){
		if (cur_well == 0){
			// Indicates that the particle was in the left well
			no_left++;
		}

		// BAOAB step
		x = BAOAB_limit(x, params, Poten_deriv, R); // Changes x according to the BAOAB limit method
		R[0] = R[1]; // Current value of R becomes the next one
		box_muller_rand(normal_dist);
		R[1] = normal_dist[0]; // Next value for R is drawn from a normal distribution

		// Recalibrate the wells if a particle has managed to cross over the barrier
		if ((cur_well == 0) && (x > 0)){
			cur_well = 1;
		}
		else if ((cur_well == 1) && (x < 0)){
			cur_well = 0;
		}


		// Attempts a lattice switch
		if (m % params->switch_regularity == 0){
			double dis = well_dis(params, cur_well, x); // Displacement from the current well
			int oth_well = (cur_well + 1) % 2; // Other well
			double diff_poten = (*Poten_shift)(x_pos(params, oth_well, dis)) - (*Poten_shift)(x); // Difference in potential
			// if uniform_rand() < min(1, exp(-diff_poten)){

			// }
		}

	}
}



int main(void){
	// double ret_arr[] = {0, 0};
	// calc_energy_diff(ret_arr, "KT", 1000000, 0.01, 1, 0.01);
	// printf("%f %f\n", ret_arr[0], ret_arr[1]);
	srand ( time(NULL) );
	char piss[] = "hello boy";
	char lips[] = " is it";
	strcat(piss, lips);
	printf("%s\n", piss);
}
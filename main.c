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
	char potential_name[10];
	long tot_timesteps;
	double timestep;
	int start_well;
	int switch_regularity;
	int write_regularity;
	
	// Physical parameters
	double kT;
	double mass;
	double shift_value;
	double minima[2];

	PotentialFun Poten;
	PotentialFun Poten_shifted;
	PotentialFun Poten_deriv;
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
void create_energy_diff_data(parameters *params, char *outputfilename){
	srand ( time(NULL) ); // Seeds the random number generators
	
	double x = x_pos(params, params->start_well, 0); // x-position initially at the bottom of the starting well

	double no_left = 0; // Number of timesteps that the particle is in the left well

	int cur_well = params->start_well; // Indicates which well the particle is in (0 is left well, 1 is right well)

	remove(outputfilename);
	FILE *outputfile;

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
		x = BAOAB_limit(x, params, params->Poten_deriv, R); // Changes x according to the BAOAB limit method
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
			double diff_poten = (*params->Poten_shifted)(x_pos(params, oth_well, dis)) - \
				(*params->Poten_shifted)(x); // Difference in potential
			// Attempts a Monte-Carlo lattice switch
			if (uniform_rand() < min(1, exp(-diff_poten))){
				cur_well = oth_well;
				x = x_pos(params, cur_well, dis);
			}

			if ((m / params->switch_regularity) % params->write_regularity == 0){
				double energy_diff = -params->kT * log((double)(no_left) / (m - no_left)) + params->shift_value;
				outputfile = fopen(outputfilename, "a");
				fprintf(outputfile, "%g\n", energy_diff);
				fclose(outputfile);
			}
		}
	}
}

/* Function to use fscanf to read strings from a line in an input file */
void read_char(FILE *input_file, char *char_value){
    /* Attempts to read line of input file */
    if (fscanf(input_file, "%s", char_value) != 1) {
        printf("Failed to read parameter\n");
        exit(1);
    }
}

/* Function to use fscanf to read doubles from a line in an input file */
void read_double(FILE *input_file, double *double_value){
    /* Attempts to read line of input file */
    if (fscanf(input_file, "%lf", double_value) != 1) {
        printf("Failed to read parameter\n");
        exit(1);
    }
}


/* Function to use fscanf to read longs from a line in an input file */
void read_long(FILE *input_file, long *long_value){
    /* Attempts to read line of input file */
    if (fscanf(input_file, "%ld", long_value) != 1) {
        printf("Failed to read parameter\n");
        exit(1);
    }
}

// Function to use fscanf to read ints from a line in an input file
void read_int(FILE *input_file, int *int_value){
	/* Attempts to read line of input file */
    if (fscanf(input_file, "%d", int_value) != 1) {
        printf("Failed to read parameter\n");
        exit(1);
    }
}

// Reads an input file to store the parameters
void store_parameters(parameters *params, char *input_filename){
	FILE *input_file = fopen(input_filename, "r");
	if (input_filename == NULL){
		printf("Failed to open input file \n");
		exit(1);
	}
	read_char(input_file,   params->potential_name);
	read_long(input_file,   &params->tot_timesteps);
	read_double(input_file, &params->timestep);
	read_int(input_file,    &params->start_well);
	read_int(input_file,    &params->switch_regularity);
	read_int(input_file,    &params->write_regularity);
	read_double(input_file, &params->kT);
	read_double(input_file, &params->mass);

	double const_arr[] = {0, 0, 0}; // Array for constants
	PotentialFun func_arr[] = {0, 0, 0}; // Array for potential functions 
	U_selector(const_arr, func_arr, params->potential_name); // Fills the constant and function arrays

	params->Poten = func_arr[0];
	params->Poten_shifted = func_arr[1];
	params->Poten_deriv = func_arr[2];

	params->minima[0] = const_arr[0]; // x-coordinates of the minima of the wells
	params->minima[1] = const_arr[1];
	params->shift_value = const_arr[2]; // The amount the right minima has been shifted upwards
}

void calculate_energy_diff(double ret_arr[2], char *output_filename, double sample_portion, int sample_regularity){
	double line_val;
	FILE *output_file = fopen(output_filename, "r");
	long line_no = 1;



}

int main(int argc, char **argv){
	char *input_filename = argv[1];
	char *output_filename = argv[2];

	parameters params;
	store_parameters(&params, input_filename);
	create_energy_diff_data(&params, output_filename);


	// long lines=0;
	// char ch;
	// FILE *fp=fopen(output_filename,"r");
	// while((ch=fgetc(fp))!=EOF)
	// {
	// if (ch=='\n') { lines++; }
	// }
	// fclose(fp);
	// printf("%ld\n", lines);


	// printf("%s\n", params.potential_name);
	// printf("%ld\n", params.tot_timesteps);
	// printf("%lf\n", params.timestep);
	// printf("%d\n", params.start_well);
	// printf("%d\n", params.switch_regularity);
	// printf("%d\n", params.write_regularity);
	// printf("%lf\n", params.kT);
	// printf("%lf\n", params.mass);
	// printf("%lf\n", params.minima[0]);
	// printf("%lf\n", params.minima[1]);
	// printf("%lf\n", params.shift_value);

	return 0;
}
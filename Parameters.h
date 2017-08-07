/*  
	Parameters.h
	Header file for dealing with the parameters which need to be passed into the lattice
	switch code. Includes a struct to store parameters, functions to read variables from
	file and a function to store parameters in a file into a "parameters" struct.
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Potentials.h"

/* Structure to store relevant parameters */
struct parameters{
	/* Numerical parameters */
	char dynamics_type[256];
	char potential_name[256];
	long tot_steps;
	int start_well;
	int switch_regularity;
	int write_regularity;
	
	/* Physical parameters */
	double kT;
	double mass;
	double shift_value;
	double minima[2];
	
	// BAOAB and BAOAB limit parameters
	double timestep;
	double R[2];

	// // Specific to the regular BAOAB method
	// double friction_param;
	// double const1;
	// double const2;
	// double const3;

	// Specific to the Monte-Carlo method
	double jump_size;

	PotentialFun Poten;
	PotentialFun Poten_shifted;
	PotentialFun Poten_deriv;
};

typedef struct parameters parameters;




/* Function to use fscanf to read strings from a line in an input file */
void read_char(FILE *input_file, char char_value[256]){
    /* Attempts to read line of input file */
    if (fscanf(input_file, "%[^\n]%*c", char_value) != 1){
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
	if (input_file == NULL){
		printf("Failed to open input file \n");
		exit(1);
	}
	read_char(input_file,   params->dynamics_type);
	read_char(input_file,   params->potential_name);
	read_long(input_file,   &params->tot_steps);
	read_int(input_file,    &params->start_well);
	read_int(input_file,    &params->switch_regularity);
	read_int(input_file,    &params->write_regularity);
	read_double(input_file, &params->kT);
	read_double(input_file, &params->mass);

	
	if (strcmp(params->dynamics_type, "BAOAB_LIMIT") == 0){
		read_double(input_file, &params->timestep);

	}
	// else if (strcmp(params->dynamics_type, "BAOAB_REGULAR") == 0){
	// 	read_double(input_file, &params->timestep);
	// 	read_double(input_file, &params->friction_param);
	// 	params->const1 = exp(-params->friction_param * params->timestep);
	// 	params->const2 = (1 - params->const1) / params->friction_param;
	// 	params->const3 = sqrt(params->kT * (1 - params->const1 * params->const1));
	// }
	else if (strcmp(params->dynamics_type, "MONTE-CARLO") == 0){
		read_double(input_file, &params->jump_size);
	}

	
	fclose(input_file);

	double poten_const_arr[] = {0, 0, 0}; // Array for potential specific constants
	PotentialFun func_arr[] = {0, 0, 0}; // Array for potential functions 
	Poten_selector(poten_const_arr, func_arr, params->potential_name); // Fills the constant and function arrays

	params->Poten = func_arr[0];
	params->Poten_shifted = func_arr[1];
	params->Poten_deriv = func_arr[2];

	params->minima[0] = poten_const_arr[0]; // x-coordinates of the minima of the wells
	params->minima[1] = poten_const_arr[1];
	params->shift_value = poten_const_arr[2]; // The amount the right minima has been shifted upwards
}

#endif
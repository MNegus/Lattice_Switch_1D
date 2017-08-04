#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "Parameters.h"
#include "Potentials.h"
#include "Dynamics.h"
#include "Random.h"
#include "MiscFunctions.h"


// Returns the x-position of a particle in space given its displacement from the current well minima
double x_pos(parameters *params, int well, double displacement){
	return displacement + params->minima[well];
}

// Returns the displacement from the current well minima given its x position in space
double well_dis(parameters *params, int well, double x_position){
	return x_position - params->minima[well];
}


double lattice_switch(double x, int cur_well, long stepno, long no_left, char *outputfilename, parameters *params){
	double dis = well_dis(params, cur_well, x); // Displacement from the current well
	int oth_well = (cur_well + 1) % 2; // Other well
	double diff_poten = (*params->Poten_shifted)(x_pos(params, oth_well, dis)) - \
		(*params->Poten_shifted)(x); // Difference in potential
	// Attempts a Monte-Carlo lattice switch
	if (uniform_rand() < min(1, exp(-diff_poten))){
		cur_well = oth_well;
		x = x_pos(params, cur_well, dis);
	}

	if ((stepno / params->switch_regularity) % params->write_regularity == 0){
		double energy_diff = -params->kT * log((double)(no_left) / (stepno - no_left)) + params->shift_value;
		FILE *outputfile = fopen(outputfilename, "a");
		fprintf(outputfile, "%g\n", energy_diff);
		fclose(outputfile);
	}
	return x;
}


// Calculates the free energy different between states in the two wells of a given potential function
void create_energy_diff_data(parameters *params, char *outputfilename){
	srand ( time(NULL) ); // Seeds the random number generators
	
	double x = x_pos(params, params->start_well, 0); // x-position initially at the bottom of the starting well

	long no_left = 0; // Number of timesteps that the particle is in the left well

	int cur_well = params->start_well; // Indicates which well the particle is in (0 is left well, 1 is right well)

	remove(outputfilename);

	DynamicsFun DynFun = Dynamics_selector(params->dynamics_type);

	if (strcmp(params->dynamics_type, "BAOAB_LIMIT") == 0){
		// Generates normally distributed values for R, which stores the current value and the value at the next timestep
		params->R[0] = box_muller_rand();
		params->R[1] = box_muller_rand();
	}
	

	// Perform lattice switching method
	for (long stepno=0; stepno < params->tot_steps; stepno++){
		if (cur_well == 0){
			// Indicates that the particle was in the left well
			no_left++;
		}

		if (isnan(x) || (x == INFINITY) || (x == -INFINITY)){
			printf("Infinite x value reached\n");
		}

		// Perform dynamics step
		x = (*DynFun)(x, params);

		// Recalibrate the wells if a particle has managed to cross over the barrier
		if ((cur_well == 0) && (x > 0)){
			cur_well = 1;
		}
		else if ((cur_well == 1) && (x < 0)){
			cur_well = 0;
		}


		// Attempts a lattice switch
		if (stepno % params->switch_regularity == 0){
			x = lattice_switch(x, cur_well, stepno, no_left, outputfilename, params);
		}
	}
}

void calculate_energy_diff(double ret_arr[2], char *output_filename, double sample_portion, int sample_regularity){
	long no_lines=0;
	char ch;
	FILE *output_file = fopen(output_filename, "r");
	while((ch=fgetc(output_file))!=EOF)
	{
		if (ch=='\n') {
		 no_lines++; 
		}
	}
	fclose(output_file);

	long start_line = (long)(no_lines * (1 - sample_portion));
	while (start_line % sample_regularity != 0){
		start_line--;
	}

	long arr_length = (long)((no_lines - start_line) / sample_regularity) + 1;
	double *data_arr = malloc(sizeof(double) * (arr_length));

	double line_val;
	double mean = 0;
	output_file = fopen(output_filename, "r");
	long line = 1;

	long arr_pos = 0;
	while (fscanf(output_file, "%lf", &line_val) == 1){
		if ((line >= start_line) && (line % sample_regularity == 0)){
			data_arr[arr_pos] = line_val;
			mean += line_val;
			arr_pos++;
		}
		line++;
	}
	fclose(output_file);

	mean = mean / arr_length;
	ret_arr[0] = mean;
	
	double std_error = 0;
	for (long j = 0; j < arr_length; j++){
		std_error += (data_arr[j] - mean) * (data_arr[j] - mean);
	}

	std_error = sqrt(std_error) / arr_length;
	ret_arr[1] = std_error;

	free(data_arr);
}

int main(int argc, char **argv){
	char *input_filename = argv[1];
	char *output_filename = argv[2];
	char *datastore_filename = argv[3];

	parameters params;
	store_parameters(&params, input_filename);

	double mean_of_simulations = 0;
	double std_error_of_simulations = 0;
	for (int i = 0; i < 10; i++){
		create_energy_diff_data(&params, output_filename);
		double ret_arr[2] = {0, 0};
		calculate_energy_diff(ret_arr, output_filename, 0.01, 100);
		mean_of_simulations += ret_arr[0];
		std_error_of_simulations += ret_arr[1] * ret_arr[1];
	}

	mean_of_simulations = mean_of_simulations / 10;
	std_error_of_simulations = sqrt(std_error_of_simulations) / 10;

	FILE *datastore_file = fopen(datastore_filename, "a");
	fprintf(datastore_file, "%s, %ld, %lf, %lf, %lf, %g, %g\n",
	 params.potential_name, params.tot_steps, params.timestep, params.kT, params.mass,
	  mean_of_simulations, std_error_of_simulations);
	fclose(datastore_file);



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
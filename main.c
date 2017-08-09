#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "Parameters.h"
#include "Dynamics.h"
#include "mt19937ar.h"

// Returns the x-position of a particle in space given its displacement from the current well minima
double x_pos(parameters *params, int well, double displacement){
	return displacement + params->minima[well];
}

// Returns the displacement from the current well minima given its x position in space
double well_dis(parameters *params, int well, double x_position){
	return x_position - params->minima[well];
}

// Attempts a lattice switch from a given position in a well
double lattice_switch(double x, int *cur_well, long stepno, long no_left, char *outputfilename, parameters *params){
	double dis = well_dis(params, *cur_well, x); // Displacement from the current well
	int oth_well = (*cur_well + 1) % 2; // Other well
	double diff_poten = (*params->Poten_shifted)(x_pos(params, oth_well, dis)) - \
		(*params->Poten_shifted)(x); // Difference in potential
	// Attempts a Monte-Carlo lattice switch
	if (genrand_real1() < min(1, exp(-diff_poten / params->kT))){
		*cur_well = oth_well;
		x = x_pos(params, *cur_well, dis);
	}

	if ((stepno / params->switch_regularity) % params->write_regularity == 0){
		double energy_diff = -params->kT * log((double)(no_left) / (stepno - no_left)) + params->shift_value;
		FILE *outputfile = fopen(outputfilename, "a");
		fprintf(outputfile, "%g\n", energy_diff);
		fclose(outputfile);
	}
	return x;
}

void add_to_bins(double x, long *bins, parameters *params){
    if (x < params->x_min){
        bins[0]++;
    }
    else if (x > params->x_max){
        bins[params->nobins - 1]++;
    }
    else {
        for (long j=1; j <= params->nobins; j++){
            if (x < params->x_min + j * params->bin_width){
                bins[j - 1]++;
                break;
            }
        }
    }
}

// Calculates the free energy different between states in the two wells of a given potential function
void create_energy_diff_data(parameters *params, char *energy_diff_filename, char *bins_filename){
	init_genrand ( time(NULL) ); // Seeds the random number generators
	
	double x = x_pos(params, params->start_well, 0); // x-position initially at the bottom of the starting well
	long no_left = 0; // Number of timesteps that the particle is in the left well
	int cur_well = params->start_well; // Indicates which well the particle is in (0 is left well, 1 is right well)
    long *bins = malloc(sizeof(long) * (params->nobins)); // Array to store number of times each bin has been visited

    for (long j = 0; j < params->nobins; j++){
        bins[j] = 0;
    }

    // Removes existing output files
    remove(energy_diff_filename);
    remove(bins_filename);

	DynamicsFun DynFun = Dynamics_selector(params->dynamics_type); // Returns the desired dynamics function

	if (strcmp(params->dynamics_type, "BAOAB_LIMIT") == 0){
		// Generates normally distributed values for R, which stores the current value and the value at the next timestep
		params->R[0] = box_muller_rand();
		params->R[1] = box_muller_rand();
	}
	

	// Perform lattice switching method
	for (long stepno=1; stepno < params->tot_steps; stepno++){
        add_to_bins(x, bins, params);

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
			x = lattice_switch(x, &cur_well, stepno, no_left, energy_diff_filename, params);
		}
	}

    FILE *bins_file = fopen(bins_filename, "w");
    for (long j=1; j <= params->nobins; j++){
        double bin_x_pos = params->x_min + j * params->bin_width;
        fprintf(bins_file, "%g, %ld\n", bin_x_pos, bins[j - 1]);
    }
    fclose(bins_file);
//    for (long j = 0; j < params->nobins; j++){
//        printf("%ld\n", bins[j]);
//    }
    free(bins);

}

// Calculates the mean and standard error of the data outputted from a lattice switching procedure
void calculate_energy_diff(double ret_arr[2], char *output_filename, double sample_portion, int sample_regularity){
	long no_lines=0; // Number of lines in file

    // Counts the number of lines in the output file
	char ch;
	FILE *output_file = fopen(output_filename, "r");
	while((ch=fgetc(output_file))!=EOF)
	{
		if (ch=='\n') {
		 no_lines++; 
		}
	}
	fclose(output_file);

    // Calculates the first line of the block of data we want to keep
	long start_line = (long)(no_lines * (1 - sample_portion));
	while (start_line % sample_regularity != 0){
        // Decreases the start line until it is divisible by the sample regularity
		start_line--;
	}

	long arr_length = (long)((double)(no_lines - start_line) / sample_regularity) + 1; // Size of array to store all the data
	double *data_arr = malloc(sizeof(double) * (arr_length));

	double line_val; // Value read in by each line
	double mean = 0; // Mean of data
	output_file = fopen(output_filename, "r");
	long line = 1; // Indicates which line of the file we are on

	long arr_pos = 0; // Position in the array we are currently reading into

    // Reads values from lines in the file to the variable line_val
	while (fscanf(output_file, "%lf", &line_val) == 1){
		if ((line >= start_line) && (line % sample_regularity == 0)){
            // Saves the data once we are at the portion of the file we want, and at the chosen regularity
			data_arr[arr_pos] = line_val;
			mean += line_val;
			arr_pos++;
		}
		line++;
	}
	fclose(output_file);

	mean = mean / arr_length;
	ret_arr[0] = mean;

    // Calculates standard error
	double std_error = 0;
	for (long j = 0; j < arr_length; j++){
		std_error += (data_arr[j] - mean) * (data_arr[j] - mean);
	}

	std_error = sqrt(std_error) / arr_length;
	ret_arr[1] = std_error;

	free(data_arr);
}

int main(int argc, char **argv){
	char *input_filename = argv[1]; // Name of the parameter input file
	char *output_filename = argv[2]; // Name of the file to output the free energy at every step
	char *datastore_filename = argv[3]; // Name of the file to store final calculated data

	parameters params; // Struct for storing parameters
	store_parameters(&params, input_filename); // Fills the parameters struct with data from input file

	double mean_of_simulations = 0;
	double std_error_of_simulations = 0;

    /* Creates file names for outputting the energies and the bins */
    char energy_diff_filename[50];
    strcpy(energy_diff_filename,  output_filename);
    strcat(energy_diff_filename, "_energies.csv");

    char bins_filename[50];
    strcpy(bins_filename, output_filename);
    strcat(bins_filename, "_bins.csv");

    // Runs the procedure 10 times to calculate an average
	for (int i = 0; i < 10; i++){
		create_energy_diff_data(&params, energy_diff_filename, bins_filename); // Runs the lattice switch procedure to create data in file
		double ret_arr[2] = {0, 0}; // Array to store average and standard error
		calculate_energy_diff(ret_arr, energy_diff_filename, 0.01, 100); // Calculates mean and standard error of the procedure
		mean_of_simulations += ret_arr[0];
		std_error_of_simulations += ret_arr[1] * ret_arr[1];
	}

	mean_of_simulations = mean_of_simulations / 10;
	std_error_of_simulations = sqrt(std_error_of_simulations) / 10;

	FILE *datastore_file = fopen(datastore_filename, "a");
    // Appends the data file with the mean and standard error, also listing the parameters associated with the run
	fprintf(datastore_file, "%s, %ld, %lf, %lf, %lf, %g, %g\n", params.potential_name, params.tot_steps,
            params.timestep, params.kT, params.mass, mean_of_simulations, std_error_of_simulations);
	fclose(datastore_file);
	return 0;
}
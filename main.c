#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "Parameters.h"
#include "Dynamics.h"
#include "mt19937ar.h"

// Returns the x-position of a particle in space given its displacement from the current well minima
double x_pos(double displacement, int well, parameters *params) {
    return displacement + params->minima[well];
}

// Returns the displacement from the current well minima given its x position in space
double well_dis(double x_position, int well, parameters *params) {
    return x_position - params->minima[well];
}

// Attempts a lattice switch from a given position in a well
double lattice_switch(double x, int *cur_well, parameters *params) {
    double dis = well_dis(x, *cur_well, params); // Displacement from the current well
    int oth_well = (*cur_well + 1) % 2; // Other well
    double diff_poten = (*params->Poten_shifted)(x_pos(dis, oth_well, params)) - \
        (*params->Poten_shifted)(x); // Difference in potential
    // Attempts a Monte-Carlo lattice switch
    if (genrand_real1() < min(1, exp(-diff_poten / params->kT))) {
        *cur_well = oth_well;
        x = x_pos(dis, *cur_well, params);
    }
    return x;
}

void add_to_bins(double x, long *bins, parameters *params) {
    if (x < params->x_min) {
        bins[0]++;
    } else if (x > params->x_max) {
        bins[params->nobins - 1]++;
    } else {
        for (long j = 1; j <= params->nobins; j++) {
            if (x < params->x_min + j * params->bin_width) {
                bins[j - 1]++;
                break;
            }
        }
    }
}



// Calculates the free energy different between states in the two wells of a given potential function
double calc_energy_difference(int savebins, char *bins_filename, parameters *params) {
    double x = x_pos(0, params->start_well, params); // x-position initially at the bottom of the starting well
    long no_left = 0; // Number of timesteps that the particle is in the left well
    int cur_well = params->start_well; // Indicates which well the particle is in (0 is left well, 1 is right well)

    long *bins; // Array for storing amount of times the walker has visited each bin

    if (savebins) {
        bins = malloc(sizeof(long) * (params->nobins)); // Array to store number of times each bin has been visited
        for (long j = 0; j < params->nobins; j++) bins[j] = 0;

        // Removes existing bin files
        remove(bins_filename);
    }

    DynamicsFun DynFun = Dynamics_selector(params->dynamics_type); // Returns the desired dynamics function

    if (strcmp(params->dynamics_type, "BAOAB_LIMIT") == 0) {
        // Generates normally distributed values for R, which stores the current value and the value at the next timestep
        params->R[0] = box_muller_rand();
        params->R[1] = box_muller_rand();
    }


    // Perform lattice switching method
    for (long stepno = 1; stepno < params->tot_steps; stepno++) {

        if (savebins) add_to_bins(x, bins, params);

        if (cur_well == 0) {
            // Indicates that the particle was in the left well
            no_left++;
        }

        if (isnan(x) || (x == INFINITY) || (x == -INFINITY)) {
            printf("Infinite x value reached\n");
        }

        // Perform dynamics step
        x = (*DynFun)(x, params);

        // Recalibrate the wells if a particle has managed to cross over the barrier
        if ((cur_well == 0) && (x > 0)) {
            cur_well = 1;
        } else if ((cur_well == 1) && (x < 0)) {
            cur_well = 0;
        }


        // Attempts a lattice switch
        if (stepno % params->switch_regularity == 0) {
            x = lattice_switch(x, &cur_well, params);
        }
    }

    if (savebins) {
        FILE *bins_file = fopen(bins_filename, "w"); // Data file to store bin data
        for (long j = 1; j <= params->nobins; j++) {
            // Fills data file with bin positions and how many hits they have
            double bin_x_pos = params->x_min + j * params->bin_width;
            fprintf(bins_file, "%g, %ld\n", bin_x_pos, bins[j - 1]);
        }
        fclose(bins_file);

        free(bins);
    }

    return -params->kT * log((double) (no_left) / (params->tot_steps - no_left)) + params->shift_value;
}

int main(int argc, char **argv) {
    char *input_filename = argv[1]; // Name of the parameter input file
    char *datastore_filename = argv[2]; // Name of the file to store final calculated data
    char *bins_filename = argv[3]; // Name of the bin output file
    char *seed_str = argv[4]; // Seed of the simulation

    char *ptr;
    unsigned long seed = (unsigned long) strtol(seed_str, &ptr, 10);

    // ///////////////////////////////
    // CHANGE WHEN ACTUALLY SIMULATING
//    seed = (unsigned long) time(NULL);
    // ///////////////////////////////

    int savebins = 1; // Variable for indicating if the data for bins will be saved or not
    if (strcmp(bins_filename, "NOBINS") == 0) savebins = 0;


    init_genrand(seed); // Seeds the random number generators

    parameters params; // Struct for storing parameters
    store_parameters(&params, input_filename); // Fills the parameters struct with data from input file

    double mean_energy_diff = 0;
    double std_error = 0;

    double *energy_differences = malloc(sizeof(double) * 10);

    // Runs the procedure 10 times to calculate an average
    for (int i = 0; i < 9; i++) {
        energy_differences[i] = calc_energy_difference(0, bins_filename, &params); // Runs the lattice switch procedure to create data in file
        mean_energy_diff += energy_differences[i];
    }

    energy_differences[9] = calc_energy_difference(savebins, bins_filename, &params); // Runs the lattice switch procedure to create data in file

    mean_energy_diff += energy_differences[9];
    mean_energy_diff /= 10;
    for (int i = 0; i < 10; i++) {
        std_error += (energy_differences[i] - mean_energy_diff) * (energy_differences[i] - mean_energy_diff);
    }
    std_error = sqrt(std_error) / 10;
    free(energy_differences);

    FILE *datastore_file = fopen(datastore_filename, "a");
    // Appends the data file with the mean and standard error, also listing the parameters associated with the run
    if (strcmp(params.dynamics_type, "BAOAB_LIMIT") == 0) {
        fprintf(datastore_file, "%s, %s, %ld, %lf, %lf, %g, %g\n", params.potential_name, params.dynamics_type,
                params.tot_steps, params.timestep, params.kT, mean_energy_diff, std_error);
    } else if (strcmp(params.dynamics_type, "MONTE-CARLO") == 0) {
        fprintf(datastore_file, "%s, %s, %ld, %lf, %lf, %g, %g\n", params.potential_name, params.dynamics_type,
                params.tot_steps, params.jump_size, params.kT, mean_energy_diff, std_error);
    }

    fclose(datastore_file);
    return 0;
}
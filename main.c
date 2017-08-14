#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
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

void calibrate_well(int *well, double x) {
    if (x < 0) *well = 0;
    else *well = 1;
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

long get_bin_no(double position, double min_bin_pos, double max_bin_pos, double bin_width, long nobins) {
    long bin_no;
    if (position < min_bin_pos) {
        bin_no = 0;
    } else if (position > max_bin_pos) {
        bin_no = nobins - 1;
    } else {
        for (long j = 1; j <= nobins; j++) {
            if (position < min_bin_pos + j * bin_width) {
                bin_no = j - 1;
                break;
            }
        }
    }

    return bin_no;
}

void add_to_bins(double x, long *bins, parameters *params) {
    long bin_no = get_bin_no(x, params->x_min, params->x_max, params->bin_width, params->nobins);
    bins[bin_no]++;
}

// The difference in the potential at a given displacement away from each well
double potential_difference(double displacement, parameters *params) {
    return (*params->Poten_shifted)(params->minima[0] + displacement) - \
     (*params->Poten_shifted)(params->minima[1] + displacement);
}

// The derivative of the difference in the potential at a given displacement away from each well
double potential_difference_deriv(double displacement, parameters *params) {
    return (*params->Poten_deriv)(params->minima[0] + displacement) - \
     (*params->Poten_deriv)(params->minima[1] + displacement);
}


// Calculates the free energy different between states in the two wells of a given potential function
double calc_free_energy_difference(int savebins, char *bins_filename, parameters *params) {
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

//
double biasing_function(double x, int cur_well, parameters *params) {
    double ret_val = 0;
    double dis = well_dis(x, cur_well, params);
    for (long gauss_num = 0; gauss_num < params->no_gaussians; gauss_num++) {
        double local_diff = potential_difference(dis, params) - params->gauss_positions[gauss_num];
        ret_val += params->height_arr[gauss_num] * exp(-params->gauss_width * local_diff * local_diff);
    }
    return ret_val;
}

double biasing_function_deriv(double x, int cur_well, parameters *params) {
    double ret_val = 0;
    double dis = well_dis(x, cur_well, params);
    for (long gauss_num = 0; gauss_num < params->no_gaussians; gauss_num++) {
        double local_diff = potential_difference(dis, params) - params->gauss_positions[gauss_num];
        ret_val += -2 * params->height_arr[gauss_num] * params->gauss_width * local_diff * \
         exp(-params->gauss_width * local_diff * local_diff) * potential_difference_deriv(dis, params);
    }
    return ret_val;
}

double biased_BAOAB_limit(double x, int cur_well, parameters *params) {
    x = x - params->timestep * \
    ((*params->Poten_deriv)(x) + params->kT * biasing_function_deriv(x, cur_well, params)) / params->mass + \
        sqrt(0.5 * params->kT * params->timestep / params->mass) * (params->R[0] + params->R[1]);

    params->R[0] = params->R[1]; // Current value of R becomes the next one
    params->R[1] = box_muller_rand(); // Next value for R is drawn from a normal distribution

    return x;

}


// Implements the lattice switch method but with biasing
double biased_simulation(double initial_f, double min_f, double gauss_width, long no_biasbins, parameters *params) {
    double x = x_pos(0, params->start_well, params); // x-position initially at the bottom of the starting well
    int cur_well = params->start_well; // Indicates which well the particle is in (0 is left well, 1 is right well)

    double f = initial_f; // f is variable used to decrease the Gaussian heights
    double cur_height = params->kT * log(f); // Height of the next Gaussian to be put down

    params->gauss_width = gauss_width;
    params->max_arrays_length = 1000; // Current length of height and Gaussian position arrays, let to be large
    params->height_arr = malloc(sizeof(double) * params->max_arrays_length); // Array to store the heights of each Gaussian
    params->gauss_positions = malloc(
            sizeof(double) * params->max_arrays_length); // Array to store the positions of each Gaussian
    params->no_gaussians = 0; // Number of Gaussians that have currently been put down

    params->pot_diff_histogram = malloc(
            sizeof(long) * no_biasbins); // Histogram for where the walker visits in potential_difference space
    for (long bin_no = 0; bin_no < no_biasbins; bin_no++) params->pot_diff_histogram[bin_no] = 0;


    /* Finds the bounds for the potential difference space we are interested in */
    double max_potential_difference = -DBL_MAX;
    double min_potential_difference = DBL_MAX;

    double min_displacement = -params->minima[1]; // Displacement can be as far left until the right well meets the origin
    double max_displacement = -params->minima[0]; // Displacement can be as far right until the left well meets the origin

    double dx = (max_displacement - min_displacement) / 100000; // Increment in x to be used in the following loop
    for (double dis = min_displacement; dis <= max_displacement; dis += dx) {
        double cur_poten_diff = potential_difference(dis,
                                                     params); // Difference in the wells given current displacements
        if (cur_poten_diff > max_potential_difference) {
            max_potential_difference = cur_poten_diff;
        }
        if (cur_poten_diff < min_potential_difference) {
            min_potential_difference = cur_poten_diff;
        }
    }
    printf("Minimum potential difference = %lf\n", min_potential_difference);

    double bin_width = (max_potential_difference - min_potential_difference) /
                       no_biasbins; // Width of bins in potential difference space


    params->R[0] = box_muller_rand();
    params->R[1] = box_muller_rand();

    long HIST_ATTEMPTS = 0;

    long gauss_regularity = 20; // Number of timesteps between each Gaussian being placed

    double absolute_minimum_x = 10;
    double absolute_maximum_x = -10;

    /* Finds the optimum set of Gaussian positions and heights to give a flat biasing distribution */
    while (f > min_f) {
        int flat_histogram = 0; // Used as Boolean to indiciate when we have a flat histogram
        long tot_steps = 0; // Total number of timesteps with current histogram
        while (flat_histogram != 1) {
            double cur_poten_diff; // Variable to store potential diff at current timestep


            HIST_ATTEMPTS++;
            if (HIST_ATTEMPTS % 1 == 0) {
                FILE *bias_out = fopen("biased_data.csv", "w");
                printf("%ld\n", HIST_ATTEMPTS);
                double x_min = -2.2;
                double x_max = -1.8;
                double dx = (x_max - x_min) / 10000;
                double cur_x = x_min;
                int cur_well;
                calibrate_well(&cur_well, cur_x);
                double biased_val = (*params->Poten_shifted)(cur_x) + biasing_function(x, cur_well, params);
                for (long i = 0; i < 10000; i++){
                    fprintf(bias_out, "%lf, %lf\n", cur_x, biased_val);
                    cur_x += dx;
                    calibrate_well(&cur_well, cur_x);
                    biased_val =  biasing_function(x, cur_well, params);
                }
                fclose(bias_out);
//
//
                exit(0);
            }

            // Performs a set number of timesteps
            for (long j = 0; j < gauss_regularity; j++) {
                // Adds to the potential difference histogram with the current position
                cur_poten_diff = potential_difference(well_dis(x, cur_well, params), params);
                params->pot_diff_histogram[get_bin_no(cur_poten_diff, min_potential_difference,
                                              max_potential_difference, bin_width, no_biasbins)]++;

                if (genrand_real1() < 0.1) {
                    // Attempts a lattice switch with probability 0.1
                    lattice_switch(x, &cur_well, params);
                } else {
                    // Else, performs a regular dynamics step
                    x = biased_BAOAB_limit(x, cur_well, params);
                }
                calibrate_well(&cur_well, x); // Recalibrates well

                if (x > absolute_maximum_x) absolute_maximum_x = x;

                if (x < absolute_minimum_x) absolute_minimum_x = x;
            }

            printf("My x value is: %lf\n", x);
            printf("My bin number is: %ld\n", get_bin_no(cur_poten_diff, min_potential_difference, max_potential_difference, bin_width, no_biasbins));
            params->height_arr[params->no_gaussians] = cur_height;
            params->gauss_positions[params->no_gaussians] = cur_poten_diff;

            params->no_gaussians++;
            if (params->no_gaussians >= params->max_arrays_length) {
                // If we've reached the limit for the arrays
                params->height_arr = longer_array(params->height_arr, params->max_arrays_length, params->max_arrays_length + 1000);
                params->gauss_positions = longer_array(params->gauss_positions, params->max_arrays_length, params->max_arrays_length + 1000);
                params->max_arrays_length += 1000;

            }

            double hist_mean = tot_steps / no_biasbins;
            int all_bins_above_threshold = 1;
            for (long bin_no = 0; bin_no < no_biasbins; bin_no++) {
                if (params->pot_diff_histogram[bin_no] <= 0.8 * hist_mean) {
                    all_bins_above_threshold = 0;
                    break;
                }
            }
            flat_histogram = all_bins_above_threshold;
        }

        f = sqrt(f);
        for (long bin_no = 0; bin_no < no_biasbins; bin_no++) params->pot_diff_histogram[bin_no] = 0;
    }


    return 0.1;
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
    seed = (unsigned long) time(NULL);
    // ///////////////////////////////

    int savebins = 1; // Variable for indicating if the data for bins will be saved or not
    if (strcmp(bins_filename, "NOBINS") == 0) savebins = 0;


    init_genrand(seed); // Seeds the random number generators

    parameters params; // Struct for storing parameters
    store_parameters(&params, input_filename); // Fills the parameters struct with data from input file


    double init_f = exp(1 / params.kT);
    double min_f = 1.1;
    double gauss_width = 0.175;
    long no_biasbins = 1000;
    biased_simulation(init_f, min_f, gauss_width, no_biasbins, &params);


//
//    double mean_energy_diff = 0;
//    double std_error = 0;
//
//    double *energy_differences = malloc(sizeof(double) * 10);
//
//    // Runs the procedure 10 times to calculate an average
//    for (int i = 0; i < 9; i++) {
//        energy_differences[i] = calc_free_energy_difference(0, bins_filename, &params); // Runs the lattice switch procedure to create data in file
//        mean_energy_diff += energy_differences[i];
//    }
//
//    energy_differences[9] = calc_free_energy_difference(savebins, bins_filename, &params); // Runs the lattice switch procedure to create data in file
//
//    mean_energy_diff += energy_differences[9];
//    mean_energy_diff /= 10;
//    for (int i = 0; i < 10; i++) {
//        std_error += (energy_differences[i] - mean_energy_diff) * (energy_differences[i] - mean_energy_diff);
//    }
//    std_error = sqrt(std_error) / 10;
//    free(energy_differences);
//
//    FILE *datastore_file = fopen(datastore_filename, "a");
//    // Appends the data file with the mean and standard error, also listing the parameters associated with the run
//    if (strcmp(params.dynamics_type, "BAOAB_LIMIT") == 0) {
//        fprintf(datastore_file, "%s, %s, %ld, %lf, %lf, %g, %g\n", params.potential_name, params.dynamics_type,
//                params.tot_steps, params.timestep, params.kT, mean_energy_diff, std_error);
//    } else if (strcmp(params.dynamics_type, "MONTE-CARLO") == 0) {
//        fprintf(datastore_file, "%s, %s, %ld, %lf, %lf, %g, %g\n", params.potential_name, params.dynamics_type,
//                params.tot_steps, params.jump_size, params.kT, mean_energy_diff, std_error);
//    }
//
//    fclose(datastore_file);
    return 0;
}
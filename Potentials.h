/*
	Potentials.h
	Header file for defining the external potentials used in lattice switching. For
	every potential, its function, its shifted function and its derivative is defined,
	and finally a function for returning pointers to relevant functions and constants
	for each potential.
*/

#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <stdio.h>
#include <math.h>
#include <string.h>


// Typedef for a function pointer
typedef double (*PotentialFun)(double);

/* Potential function from Kinetic Theory notes */

// External potential
double KT_Poten(double x) {
    double pot_val; // Value of the potential at x
    if (x <= -1) {
        pot_val = 5 * (x + 2) * (x + 2);
    } else if (x <= 1.2) {
        pot_val = 10 - 5 * x * x;
    } else {
        pot_val = 5 * (x - 2.4) * (x - 2.4) - 4.4;
    }
    return pot_val;
}

// External potential shifted
double KT_Poten_shifted(double x) {
    double pot_val; // Value of the potential at x
    if (x <= -1) {
        pot_val = 5 * (x + 2) * (x + 2);
    } else if (x <= 0) {
        pot_val = 10 - 5 * x * x;
    } else if (x <= 1.2) {
        pot_val = 14.4 - 5 * x * x;
    } else {
        pot_val = 5 * (x - 2.4) * (x - 2.4);
    }
    return pot_val;
}

// Derivative of potential function
double KT_Poten_deriv(double x) {
    double pot_val; // Value of the potential at x
    if (x <= -1) {
        pot_val = 10 * (x + 2);
    } else if (x <= 1.2) {
        pot_val = -10 * x;
    } else {
        pot_val = 10 * (x - 2.4);
    }
    return pot_val;
}


/* A quartic potential function */

// External potential
double QUARTIC_Poten(double x) {
    double x_0 = -0.126000192586256;
    return pow(x + x_0, 4) - 4 * (x + x_0) * (x + x_0) - (x + x_0) + 2.618555980765;
}

// External potential shifted
double QUARTIC_Poten_shifted(double x) {
    double pot_val; // Value of the potential at x
    double x_0 = -0.126000192586256;
    if (x < 0) {
        pot_val = pow(x + x_0, 4) - 4 * (x + x_0) * (x + x_0) - (x + x_0) + 2.618555980765;
    } else {
        pot_val = pow(x + x_0, 4) - 4 * (x + x_0) * (x + x_0) - (x + x_0) + 5.444192066610897;
    }
    return pot_val;
}

// Derivative of potential function
double QUARTIC_Poten_deriv(double x) {
    double x_0 = -0.126000192586256;
    return 4 * pow(x + x_0, 3) - 8 * (x + x_0) - 1;
}


/* Potential with differing well widths */

// External potential
double DIFF_WIDTH_Poten(double x) {
    double pot_val; // Value of the potential at x
    if (x <= -1) {
        pot_val = 5 * (x + 2) * (x + 2);
    } else if (x <= 0) {
        pot_val = 10 - 5 * x * x;
    } else if (x <= 0.3461) {
        pot_val = -50 * x * x + 10;
    } else {
        pot_val = 50 * (x - sqrt(0.48)) * (x - sqrt(0.48)) - 2;
    }
    return pot_val;
}

// External potential shifted
double DIFF_WIDTH_Poten_shifted(double x) {
    double pot_val; // Value of the potential at x
    if (x <= -1) {
        pot_val = 5 * (x + 2) * (x + 2);
    } else if (x <= 0) {
        pot_val = 10 - 5 * x * x;
    } else if (x <= 0.3461) {
        pot_val = -50 * x * x + 12;
    } else {
        pot_val = 50 * (x - sqrt(0.48)) * (x - sqrt(0.48));
    }
    return pot_val;
}

// Derivative of potential function
double DIFF_WIDTH_Poten_deriv(double x) {
    double pot_val; // Value of the potential at x
    if (x <= -1) {
        pot_val = 10 * (x + 2);
    } else if (x <= 0) {
        pot_val = -10 * x;
    } else if (x <= 0.3461) {
        pot_val = -100 * x;
    } else {
        pot_val = 100 * (x - sqrt(0.48));
    }
    return pot_val;
}

// Returns potential specific constant values and function pointers
void Poten_selector(double const_arr[], PotentialFun func_arr[], char name[]) {
    // Function from kinetic theory notes
    if (strcmp(name, "KT") == 0) {
        const_arr[0] = -2; // Left minimum
        const_arr[1] = 2.4; // Right minimum
        const_arr[2] = 4.4; // Amount the right well is shifted

        // Function pointers
        func_arr[0] = &KT_Poten;
        func_arr[1] = &KT_Poten_shifted;
        func_arr[2] = &KT_Poten_deriv;
    }
        // Quartic function
    else if (strcmp(name, "QUARTIC") == 0) {
        const_arr[0] = -1.220997215942; // Left minimum
        const_arr[1] = 1.5989977937; // Right minimum
        const_arr[2] = 2.8256360858458973; // Amount the right well is shifted

        // Function pointers
        func_arr[0] = &QUARTIC_Poten;
        func_arr[1] = &QUARTIC_Poten_shifted;
        func_arr[2] = &QUARTIC_Poten_deriv;
    }
        // Potential function with differing well widths
    else if (strcmp(name, "DIFF_WIDTH") == 0) {
        const_arr[0] = -2; // Left minimum
        const_arr[1] = sqrt(0.48); // Right minimum
        const_arr[2] = 2; // Amount the right well is shifted

        // Function pointers
        func_arr[0] = &DIFF_WIDTH_Poten;
        func_arr[1] = &DIFF_WIDTH_Poten_shifted;
        func_arr[2] = &DIFF_WIDTH_Poten_deriv;
    }
}


#endif // POTENTIALS_H
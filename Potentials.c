#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Potentials.h"

// Returns potential specific constant values and function pointers
void U_selector(double const_arr[], PotentialFun func_arr[], char name[]) {
	// Function from kinetic theory notes
	if (strcmp(name, "KT") == 0) {
		const_arr[0] = -2; // Left minimum
		const_arr[1] = 2.4; // Right minimum
		const_arr[2] = 4.4; // Amount the right well is shifted

		// Function pointers
		func_arr[0] = &KT_U;
		func_arr[1] = &KT_U_shifted;
		func_arr[2] = &KT_DU;
	} 
	// Quartic function
	else if (strcmp(name, "QUARTIC") == 0) {
		const_arr[0] = -1.220997215942; // Left minimum
		const_arr[1] = 1.5989977937; // Right minimum
		const_arr[2] = 2.8256360858458973; // Amount the right well is shifted

		// Function pointers
		func_arr[0] = &QUARTIC_U;
		func_arr[1] = &QUARTIC_U_shifted;
		func_arr[2] = &QUARTIC_DU;
	}
	// Potential function with differing well widths
	else if (strcmp(name, "DIFF_WIDTH") == 0) {
		const_arr[0] = -2; // Left minimum
		const_arr[1] = sqrt(0.48); // Right minimum
		const_arr[2] = 2; // Amount the right well is shifted

		// Function pointers
		func_arr[0] = &DIFF_WIDTH_U;
		func_arr[1] = &DIFF_WIDTH_U_shifted;
		func_arr[2] = &DIFF_WIDTH_DU;
	}
}


/* Potential function from Kinetic Theory notes */

// External potential
double KT_U(double x){
	if (x <= -1) {
		return 5 * (x + 2) * (x + 2);
	} 
	else if (x <= 1.2) {
		return 10 - 5 * x * x;
	}
	else {
		return 5 * (x - 2.4) * (x - 2.4) - 4.4;
	}
}

// External potential shifted
double KT_U_shifted(double x){
	if (x <= -1){
		return 5 * (x + 2) * (x + 2);
	}
	else if (x <= 0) {
		return 10 - 5 * x * x;
	}
	else if (x <= 1.2) {
		return 14.4 - 5 * x * x;
	}
	else {
		return 5 * (x - 2.4) * (x - 2.4);
	}
}

// Derivative of potential function
double KT_DU(double x){
	if (x <= -1) {
		return 10 * (x + 2);
	}
	else if (x <= 1.2) {
		return -10 * x;
	}
	else {
		return 10 * (x - 2.4);
	}
}


/* A quartic potential function */

// External potential
double QUARTIC_U(double x){
	double x_0 = -0.126000192586256;
	return pow(x + x_0, 4) - 4 * (x + x_0) * (x + x_0) - (x + x_0) + 2.618555980765;
}

// External potential shifted
double QUARTIC_U_shifted(double x){
	double x_0 = -0.126000192586256;
	if (x < 0) {
		return pow(x + x_0, 4) - 4 * (x + x_0) * (x + x_0) - (x + x_0) + 2.618555980765;
	}
	else {
		pow(x + x_0, 4) - 4 * (x + x_0) * (x + x_0) - (x + x_0) + 5.444192066610897;
	}
}

// Derivative of potential function
double QUARTIC_DU(double x){
	double x_0 = -0.126000192586256;
	return  4 * pow(x + x_0, 3) - 8 * (x + x_0) - 1;
}


/* Potential with differing well widths */

// External potential
double DIFF_WIDTH_U(double x){
	if (x <= -1) {
		return 5 * (x + 2) * (x + 2);
	}
	else if (x <= 0) {
		return 10 - 5 * x * x;
	}
	else if (x <= 0.3461) {
		return -50 * x * x + 10;
	}
	else {
		return 50 * (x - sqrt(0.48)) * (x - sqrt(0.48)) - 2;
	}
}

// External potential shifted
double DIFF_WIDTH_U_shifted(double x){
	if (x <= -1) {
		return 5 * (x + 2) * (x + 2);
	}
	else if (x <= 0) {
		return 10 - 5 * x * x;
	}
	else if (x <= 0.3461) {
		return -50 * x * x + 12;
	}
	else {
		return 50 * (x - sqrt(0.48)) * (x - sqrt(0.48));
	}
}

// Derivative of potential function
double DIFF_WIDTH_DU(double x){
	if (x <= -1) {
		return 10 * (x + 2);
	}
	else if (x <= 0) {
		return -10 * x;
	}
	else if (x <= 0.3461) {
		return -100 * x;
	}
	else {
		return 100 * (x - sqrt(0.48));
	}
}
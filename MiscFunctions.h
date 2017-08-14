/*
	MiscFunctions.h
	Header file for miscellaneous functions
*/

#ifndef MISCFUNCTIONS_H
#define MISCFUNCTIONS_H

#include <math.h>

// Returns the minimum of two double values
double min(double x1, double x2) {
    double ret_var; // Variable to return
    if (x1 < x2) {
        ret_var = x1;
    } else {
        ret_var = x2;
    }
    return ret_var;
}

// "Extends the length of an array" by creating a new one
double *longer_array(double *cur_arr, long cur_arr_length, long new_arr_length){
    double *temp = malloc(sizeof(double) * new_arr_length);
    for (long arr_index=0; arr_index < cur_arr_length; arr_index++){
        temp[arr_index] = cur_arr[arr_index];
    }
    return temp;
}


// Gaussian function with given width, height and mean
double gaussian(double position, double width, double mean, double height){
    return height * exp(-width * (position - mean) * (position - mean));
}

// Derivative of a gaussian with given width, height and mean
double gaussian_deriv(double position, double width, double mean, double height){
    return -2 * width * (position - mean) * gaussian(position, width, mean, height);
}

#endif
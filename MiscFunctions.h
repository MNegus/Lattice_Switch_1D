/*
	MiscFunctions.h
	Header file for miscellaneous functions
*/

#ifndef MISCFUNCTIONS_H
#define MISCFUNCTIONS_H

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

#endif
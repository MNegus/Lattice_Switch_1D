/*
	MiscFunctions.h
	Header file for miscellaneous functions
*/

#ifndef MISCFUNCTIONS_H
#define MISCFUNCTIONS_H

// Returns the minimum of two double values
double min(double x1, double x2){
	if (x1 < x2) {
		return x1;
	}
	else {
		return x2;
	}
}

#endif
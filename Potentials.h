#ifndef POTENTIALS_H
#define POTENTIALS_H

// Typedef for a function pointer
typedef double (*PotentialFun)(double);

// Returns potential specific constant values and function pointers
void U_selector(double const_arr[], PotentialFun func_arr[], char name[]);

// Potential function from Kinetic Theory notes
double KT_U(double x);
double KT_U_shifted(double x);
double KT_DU(double x);

// Quartic potential
double QUARTIC_U(double x);
double QUARTIC_U_shifted(double x);
double QUARTIC_DU(double x);

// Potential with differing well widths
double DIFF_WIDTH_U(double x);
double DIFF_WIDTH_U_shifted(double x);
double DIFF_WIDTH_DU(double x);


#endif // POTENTIALS_H
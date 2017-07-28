#include <stdio.h>
#include "Potentials.h"


int main(void){
	double const_arr[] = {0, 0, 0};
	PotentialDef func_arr[] = {0, 0, 0};
	U_selector(const_arr, func_arr, "DIFF_WIDTH");
	printf("%f %f %f\n", const_arr[0], const_arr[1], const_arr[2]);
	return 0;
}
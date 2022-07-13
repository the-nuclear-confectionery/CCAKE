#include <stdio.h>

#include "lib.h"

int main(void)
{
	printf("Starting!\n");
	initialize("../Coefficients_Parameters.dat");

	double point[4] = {200.0, 50.0, 75.0, 100.0};
	double densities[4];
	get_densities(point, densities);
	printf("Found these densities: %lf %lf %lf %lf\n",
		densities[0], densities[1], densities[2], densities[3]);

	printf("Finished!\n");

	return 0;
}

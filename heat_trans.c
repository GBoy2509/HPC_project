// finite-difference method to slove 2D heat transfer problem
// Square Omega of dimension 1 by 1
// Euler scheme for the time discretization
// finite difference method for the spatial discretization
//
#include <stdio.h>
#include "global.h"

int main()
{
	float dt = Tend / (float)nsteps;
	float dx = (float)L / (float)N;
	float dy = (float)L / (float)N;
	printf ("%d, %f, %f, %f, %d\n", N, dx, dy, dt, nsteps);

	return 0;
}

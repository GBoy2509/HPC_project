// finite-difference method to slove 2D heat transfer problem
// Square Omega of dimension 1 by 1
// Euler scheme for the time discretization
// finite difference method for the spatial discretization
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global.h"

const double pi = 4.0 * atan(1.0);
int main() {
	double u[N+1][N+1], unew[N+1][N+1];
	double x, y;
	int it;
	double dt;
	dt = Tend / (double)nsteps;
	double dx;
	dx = (double)L / (double)N;
	double dy;
	dy = (double)L / (double)N;
	int i,j,k,iter;
	FILE *init, *u01;
	
	double forward(double ujk, double ujmk, double ujpk, double ujkm, double ujkp, double x, double y);
	printf ("%d, %f, %f, %f, %d\n", N, dx, dy, dt, nsteps);
	
	// grid init (x-d first, then y-d)
	// use initial condition
	for (j = 1; j < (N+1); j++) {
		y = (j - 0.5)*dy;
		for (i = 1; i < (N+1); i++) {
			x = (i - 0.5)*dx;
			u[i][j] = exp(x+y);
		}
	}
	// use boundary conditions
	// just set top & bottom as tau_g, left & right as tau_h
	for (i = 1; i < (N+1); i++) {
		u[i][0] = g;
		u[i][N+1] = g;
	}
	for (i = 1; i < (N+1); i++) {
		u[N+1][i] = u[N][i] + dx*h/k;
		u[0][i] = u[1][i] + dx*h/k;
	}
	// test code
	init = fopen("init.csv","w+");
	for (j = 0; j < (N+2); j++) {
		y = (j - 0.5)*dy;
                for (i = 0; i < (N+2); i++) {
			x = (i - 0.5)*dx;
                        fwrite(init,"%f,%f,%f\n",&x,&y,&u[i][j]);
                }
        }
	fclose(init);

	// calculate next step
	// forward euler scheme
	for (iter = 1; iter < (nsteps+1); iter++) {
		for (k = 1; k < (N+1); k++) {
			y = (k - 0.5)*dy;
        	        for (j = 1; j < (N+1); j++) {
				x = (j - 0.5)*dx;
				unew[j][k] = forward(u[j][k],u[j-1][k],u[j+1][k],u[j][k-1],u[j][k+1],x,y);
			}
		}
	}
	for (i = 1; i < (N+1); i++) {
                unew[i][0] = g;
                unew[i][N+1] = g;
        }
        for (i = 1; i < (N+1); i++) {
                unew[N+1][i] = unew[N][i] + dx*h/k;
                unew[0][i] = unew[1][i] + dx*h/k;
        }
	// backward euler is constructing
	//
	//
	// update
	for (j = 0; j < (N+2); j++) {
                for (i = 0; i < (N+2); i++) {
			u[i][j] = unew[i][j];
                }
        }

	// output
	// change as your needs (here, I print out t~=0.1)
	// demo ver with csv file, should change with HDF5!!!
	if (it == (int)(0.1/dt)) {
		u01 = fopen("u01.csv","w+");
        		for (j = 0; j < (N+2); j++) {
                		y = (j - 0.5)*dy;
                			for (i = 0; i < (N+2); i++) {
                        			x = (i - 0.5)*dx;
                        			fwrite(u01,"%f,%f,%f\n",&x,&y,&u[i][j]);
                			}
        		}
	}
        fclose(u01);
	return 0;
}

// euler scheme
double forward(double ujk, double ujmk, double ujpk, double ujkm, double ujkp, double x, double y) {
        double u1, u2;
	double ujknew;
        double dt;
	double f;
	int i,j;

        dt = Tend / (double)nsteps;
        double dx;
        dx = (double)L / (double)N;
        double dy;
        dy = (double)L / (double)N;
	
        double CFL1 = k*dt/(dx*dx);
        double CFL2 = k*dt/(dy*dy);
	f = sin(L*pi*(x+y));
        u1 = CFL1*(ujpk-2*ujk+ujmk);
        u2 = CFL2*(ujkp-2*ujk+ujkm);
        ujknew = ujk + 1/(rho*c)*(f+u1+u2);

        return ujknew;
}



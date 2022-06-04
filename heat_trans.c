// finite-difference method to slove 2D heat transfer problem
// Square Omega of dimension 1 by 1
// Euler scheme for the time discretization
// finite difference method for the spatial discretization
// use HDF5 to store output data
//
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "global.h"
#include "hdf5.h"
#define FILE2 "heat_trans_data.h5"


int main() {

	// allocate two dimensional array (N+1)*(N+1)
	double** u = (double**)malloc(sizeof(float*) * (N + 2));
	for (int i = 0; i < (N + 2); i++) {
		u[i] = (double*)malloc(sizeof(double) * (N + 2));
	}
	double** unew = (double**)malloc(sizeof(float*) * (N + 2));
	for (int i = 0; i < (N + 2); i++) {
		unew[i] = (double*)malloc(sizeof(double) * (N + 2));
	}
	
	double x, y;
	double dt;
	dt = Tend / (double)nsteps;
	double dx;
	dx = (double)L / (double)N;
	double dy;
	dy = (double)L / (double)N;
	FILE *init;
	FILE *uout;
	double pi = 4.0 * atan(1.0);

	printf ("%d, %f, %f, %f, %d\n", N, dx, dy, dt, nsteps);
	
	double forward(double ujk, double ujmk, double ujpk, double ujkm, double ujkp, double x, double y);

	u[0][0] = 0.0;
	unew[0][0] = 0.0;
	u[N+1][N+1] = 0.0;
	unew[N+1][N+1] = 0.0;
	u[0][N+1] = 0.0;
	unew[0][N+1] = 0.0;
	u[N+1][0] = 0.0;
	unew[N+1][0] = 0.0;
	// grid init (x-d first, then y-d)
	// use initial condition
	for (int j = 1; j < (N+1); j++) {
		y = (j - 0.5)*dy;
		for (int i = 1; i < (N+1); i++) {
			x = (i - 0.5)*dx;
			u[i][j] = exp(x+y);
		}
	}
	// use boundary conditions
	// just set top & bottom as tau_g, left & right as tau_h
	for (int i = 1; i < (N+1); i++) {
		u[i][0] = g;
		u[i][N+1] = g;
	}
	for (int i = 1; i < (N+1); i++) {
		u[N+1][i] = u[N][i] + dx*h/kcond;
		u[0][i] = u[1][i] + dx*h/kcond;
	}
	// test code
	init = fopen("init.csv","w+");
	for (int j = 0; j < (N+2); j++) {
		y = (j - 0.5)*dy;
        	for (int i = 0; i < (N+2); i++) {
			x = (i - 0.5)*dx;
            		fprintf(init,"%f,%f,%f\n",x,y,u[i][j]);
                }
    	}
	fclose(init);

	// calculate next step
	// forward euler scheme
	for (int iter = 1; iter < (nsteps + 1); iter++) {
		for (int k = 1; k < (N + 1); k++) {
			y = (k - 0.5) * dy;
			for (int j = 1; j < (N + 1); j++) {
				x = (j - 0.5) * dx;
				unew[j][k] = forward(u[j][k], u[j - 1][k], u[j + 1][k], u[j][k - 1], u[j][k + 1], x, y);
			}
		}

		for (int i = 1; i < (N + 1); i++) {
			unew[i][0] = g;
			unew[i][N + 1] = g;
		}
		for (int i = 1; i < (N + 1); i++) {
			unew[N + 1][i] = unew[N][i] + dx * h / kcond;
			unew[0][i] = unew[1][i] + dx * h / kcond;
		}
		// backward euler is constructing
		//
		//
		// update
		for (int j = 0; j < (N + 2); j++) {
			for (int i = 0; i < (N + 2); i++) {
				u[i][j] = unew[i][j];
			}
		}

		// output
		// change as your needs (here, I print out t~=0.01)
		// demo ver with csv file, should change with HDF5!!!
		// still have problems here, should u get higher constantly in this BCs ?!!!
		if (iter == (int)(0.01 / dt)) {
			uout = fopen("uout.csv", "w+");
			for (int j = 0; j < (N + 2); j++) {
				y = (j - 0.5) * dy;
				for (int i = 0; i < (N + 2); i++) {
					x = (i - 0.5) * dx;
					fprintf(uout, "%f,%f,%f\n", x, y, u[i][j]);
				}
			}
			fclose(uout);
		}
		

		iter = iter + 1;
	}
	for (int i = 0; i < (N + 2); i++) {
		free(u[i]);
		free(unew[i]);
	}
	free(u);
	free(unew);



    /* HDF5 initialization */

    hid_t        file_id, dataset_id, group_id, dataspace_id;  /* identifiers */
    hsize_t      dims[2];
    herr_t       status;
    int * vec1 = (int*)malloc((N+2)*(N+2)*sizeof(int));   
    free(vec1);

    /* HDF5: Create a new file to store velocity datasets. */
    file_id = H5Fcreate("heat_trans_data.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = N+2; 
    dims[1] = N+2;

    dataspace_id = H5Screate_simple(2, dims, NULL);

    /* Create two groups to store initial values and output values in the file. */
    group_id = H5Gcreate2(file_id, "/output", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Create the datasets. */
    dataset_id = H5Dcreate2(file_id, "/output/uout", H5T_STD_I32BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    //  printf("original dset_data[0][0][0]:%2d\n", dset_data[0][0][0]);

     /* Write the first dataset. */
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       u);

     /* Close the data space for the first dataset. */
    status = H5Sclose(dataspace_id);

     /* Close the first dataset. */
    status = H5Dclose(dataset_id);
     /* Close the group. */
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);


	return 0;
}


    // eulersch function
    double forward(double ujk, double ujmk, double ujpk, double ujkm, double ujkp, double x, double y) {
        double u1, u2;
        double ujknew;
        double dt;
        double f;

        dt = Tend / (double)nsteps;
        double dx;
        dx = (double)L / (double)N;
        double dy;
        dy = (double)L / (double)N;
        double pi = 4.0 * atan(1.0);

        double CFL1 = kcond*dt/(dx*dx);
        double CFL2 = kcond*dt/(dy*dy);
        f = sin(L*pi*(x+y));
        u1 = CFL1*(ujpk-2*ujk+ujmk);
        u2 = CFL2*(ujkp-2*ujk+ujkm);
        ujknew = ujk + 1/(rho*c)*(f+u1+u2);

        return ujknew;

}



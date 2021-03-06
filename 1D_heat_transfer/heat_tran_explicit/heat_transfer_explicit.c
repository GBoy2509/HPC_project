// Finite-difference method to slove 1D heat transfer problem
// Square of dimension 1 by 1
// Explicit Euler scheme for the time discretization
// Finite difference method for the spatial discretization
// Use HDF5 to store output data
// 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <petscksp.h>
#include <petsc.h>
#include <petscvec.h>
#include "petscmat.h"
#include "global.h"
#include "hdf5.h"
#define FILE2 "heat_trans_data_implicit.h5"

static char help[] = "Explicit EULER method\n\n";

int main(int argc, char** argv) {
	// allocate one dimensional array
	//double* u = (double**)malloc(sizeof(float*) * (N + 2));
	//for (int i = 0; i < (N + 2); i++) {
	//	u[i] = (double*)malloc(sizeof(double) * (N + 2));
	//}
	//double** unew = (double**)malloc(sizeof(float*) * (N + 2));
	//for (int i = 0; i < (N + 2); i++) {
	//	unew[i] = (double*)malloc(sizeof(double) * (N + 2));
	//} 
	//double** f = (double**)malloc(sizeof(float*) * (N));
	//for (int i = 0; i < (N); i++) {
	//	f[i] = (double*)malloc(sizeof(double) * (N));
	//}
	double uexact[N];
	double x;
	double dx;
	int ii;
	dx = (double)L / (double)N;
	//double dy;
	//dy = (double)L / (double)N;
	// FILE* init;
	// FILE* uout;

	// calculate analytical solution
	double pi = 4.0 * atan(1.0);
	double coef = 1.0 / kcond / L / L / pi / pi;
	double CFL = dt * kcond / dx / dx / rho / c;
	for (int i = 0; i < (N + 1); i++) {
		x = (i + 0.5) * dx;
		uexact[i] = coef * sin(L * pi * x) - coef * sin(L * pi) * x;
		printf("%f\n", uexact[i]);
	}

	// HDF5 part
	// Initialization
	hid_t        file_id, dataset_id1, dataset_id2, group_id, dataspace_id1, dataspace_id2;  /* identifiers */
	hsize_t      dim1[1], dim2[1];
	herr_t       status;
	double* vec1 = (double*)malloc(N * sizeof(double));
	double* vec2 = (double*)malloc(3 * sizeof(double));
	free(vec1);
	free(vec2);

	dim1[0] = N;
	dim2[0] = 3;
	// grid init
	// use initial condition
	//for (int i = 1; i < (N + 1); i++) {
	//	x = (i - 0.5) * dx;
	//	u[i] = exp(x);
	//}

	// use boundary conditions
	//u[0] = g;
	//u[N + 1] = g;
	
	// test code
	//init = fopen("init.csv", "w+");
	//for (int i = 0; i < (N + 2); i++) {
	//	x = (i - 0.5) * dx;
	//	fprintf(init, "%f, %f\n", x, u[i]);
	//}
	//fclose(init);

	// Petsc Part
	//MPI_Comm       comm;
	//PetscMPIInt    rank;

	Mat            A;
	Vec            u, uold, f;
	PetscInt       nlocal, rstart, rend;
	PetscInt       n = N;
	PetscInt	   col[3], num[N];
	PetscReal      temp, t = 0.0, tend = 1.0;
	//PetscReal      norm0 = 0.0, norm1 = 1.0, tor = 1.e-8, err;
	PetscScalar    value[3], uout[N], para[3];
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;

	// MPI initial
	//comm = MPI_COMM_WORLD;
	//MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);

	// Create temp vector u, uold
	ierr = VecCreate(PETSC_COMM_WORLD, &u); CHKERRQ(ierr);
	ierr = VecSetSizes(u, PETSC_DECIDE, n); CHKERRQ(ierr);
	ierr = VecSetFromOptions(u); CHKERRQ(ierr);
	/* Identify the starting and ending mesh points on each
    processor for the interior part of the mesh. We let PETSc decide
    above. */
	ierr = VecGetOwnershipRange(u, &rstart, &rend); CHKERRQ(ierr);
	ierr = VecGetLocalSize(u, &nlocal); CHKERRQ(ierr);
	ierr = VecDuplicate(u, &uold); CHKERRQ(ierr);
	ierr = VecDuplicate(u, &f); CHKERRQ(ierr);
	/*
    Create matrix.  When using MatCreate(), the matrix format can
    be specified at runtime.

    Performance tuning note:  For problems of substantial size,
    preallocation of matrix memory is crucial for attaining good
    performance. See the matrix chapter of the users manual for details.

    We pass in nlocal as the "local" size of the matrix to force it
    to have the same parallel layout as the vector created above.
	*/
	// Create coefs matrix A
	ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
	ierr = MatSetSizes(A, nlocal, nlocal, n, n); CHKERRQ(ierr);
	ierr = MatSetFromOptions(A); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(A, 3, NULL, 2, NULL); CHKERRQ(ierr);
	ierr = MatSetUp(A); CHKERRQ(ierr);

	// Set values to A
	if (!rstart)
	{
		rstart = 1;
		ii = 0; col[0] = 0; col[1] = 1; value[0] = 1.0 - 2.0 * CFL; value[1] = CFL;
		ierr = MatSetValues(A, 1, &ii, 2, col, value, INSERT_VALUES); CHKERRQ(ierr);
	}
	if (rend == n)
	{
		rend = n - 1;
		ii = n - 1; col[0] = n - 2; col[1] = n - 1; value[0] = CFL; value[1] = 1.0 - 2.0 * CFL;
		ierr = MatSetValues(A, 1, &ii, 2, col, value, INSERT_VALUES); CHKERRQ(ierr);
	}	
	value[0] = CFL;
	value[1] = 1.0 - 2.0 * CFL;
	value[2] = CFL;
	for (int i = rstart; i < rend; i++) {
		col[0] = i - 1; col[1] = i; col[2] = i + 1;
		MatSetValues(A, 1, &i, 3, col, value, INSERT_VALUES);
	}
	/* Assemble the matrix */
	ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	//ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

	// Set values to f
	ierr = VecGetOwnershipRange(f, &rstart, &rend); CHKERRQ(ierr);
	for (int i = rstart; i < rend; i++)
	{
		temp = dt / rho / c * sin(L * pi * ((double)i - 0.5) * dx);
		ierr = VecSetValue(f, i, temp, INSERT_VALUES);
	}
	ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(f); CHKERRQ(ierr);
	//ierr = VecView(f, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

	// restart=true: restart
	bool restart = false;
	if (restart == false)
	{
		/* HDF5: Create a new file to store datasets. */
		file_id = H5Fcreate(FILE2, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		/* HDF5: data1-output velocity
				 data2-parameters:t, dt, N */
		dataspace_id1 = H5Screate_simple(1, dim1, NULL);
		dataspace_id2 = H5Screate_simple(1, dim2, NULL);

		/* Create two groups to store initial values and output values in the file. */
		group_id = H5Gcreate2(file_id, "/output", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		/* Create the datasets. */
		dataset_id1 = H5Dcreate2(file_id, "/output/uout", H5T_IEEE_F64BE, dataspace_id1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		dataset_id2 = H5Dcreate2(file_id, "/output/parameters", H5T_IEEE_F64BE, dataspace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		para[0] = t;
		para[1] = dt;
		para[2] = n;

		// Set values to u
		ierr = VecGetOwnershipRange(u, &rstart, &rend); CHKERRQ(ierr);
		for (int i = rstart; i < rend; i++) {
			temp = exp(((double)i - 0.5) * dx);
			ierr = VecSetValue(u, i, temp, INSERT_VALUES); CHKERRQ(ierr);
		}
		ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
		//ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}

	if (restart == true) {
		/* HDF5: Open an existing file which stores datasets. */
		file_id = H5Fopen(FILE2, H5F_ACC_RDWR, H5P_DEFAULT);

		/* HDF5: data1-output velocity
				 data2-parameters:t, dt, N */

				 /* Create two groups to store initial values and output values in the file. */
		group_id = H5Gopen(file_id, "/output", H5P_DEFAULT);

		/* Create the datasets. */
		dataset_id1 = H5Dopen(file_id, "/output/uout", H5P_DEFAULT);
		dataset_id2 = H5Dopen(file_id, "/output/parameters", H5P_DEFAULT);

		/* Write data from hdf5 to arrays */
		status = H5Dread(dataset_id1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uout);
		status = H5Dread(dataset_id2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, para);

		t = para[0];
		dt = para[1];
		n = para[2];

		// Set values to u-from uout in hdf5
		ierr = VecGetOwnershipRange(u, &rstart, &rend); CHKERRQ(ierr);
		for (int i = rstart; i < rend; i++) {
			temp = uout[i];
			ierr = VecSetValue(u, i, temp, INSERT_VALUES); CHKERRQ(ierr);
		}
		ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
		ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
		//ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}
	// calculate next step
	// forward euler scheme
	for (int iter = 0; iter < nstep; iter++) {
		ierr = VecCopy(u, uold); CHKERRQ(ierr);
		ierr = MatMultAdd(A, uold, f, u); CHKERRQ(ierr);
		t = t + dt;
		para[0] = t;

		//ierr = VecNorm(u, NORM_1, &norm1);
		//err = fabs(norm1 - norm0);
		//if (err < tor) {
		//	PetscPrintf(comm, "Number of iteration is %d, err is %f\n", iter + 1, err);
		//	break;
		//}
		//norm0 = norm1;
		if ((iter + 1) % 10 == 0) {
			for (int i = 0; i < N; i++) {
				num[i] = i;
			}
			ierr = VecGetValues(u, N, num, uout); CHKERRQ(ierr);
			status = H5Dwrite(dataset_id1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uout);
			status = H5Dwrite(dataset_id2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, para);
		}
		if (t > tend) {
			break;
		}
		//ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}
	ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	
	// output
	// change as your needs (here, I print out t~=0.01)
	// demo ver with csv file, should change with HDF5!!!
	// still have problems here, should u get higher constantly in this BCs ?!!!
	//	if (iter == (int)(0.01 / dt)) {
	//		uout = fopen("uout.csv", "w+");
	//		for (int j = 0; j < (N + 2); j++) {
	//			y = (j - 0.5) * dy;
	//			for (int i = 0; i < (N + 2); i++) {
	//				x = (i - 0.5) * dx;
	//				fprintf(uout, "%f,%f,%f\n", x, y, u[i][j]);
	//			}
	//		}
	//		fclose(uout);
	//	}

	//	iter = iter + 1;
	//}

	for (int i = 0; i < N; i++)
	{
		num[i] = i;
	}
	ierr = VecGetValues(u, N, num, uout); CHKERRQ(ierr);
	//ierr = PetscScalarView(n, uout, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

	// calculate max error
	double err1 = 0.0, err2;
	for (int i = 0; i < N; i++) {
		err2 = fabs(uexact[i] - uout[i]);
		//printf("%lf,%lf,%lf\n", uexact[i], uout[i], err2);
		if (err2 > err1) {
			err1 = err2;
		}
		//printf("The max error is : %lf\n", err1);
	}
	printf("The final error is : %lf\n", err1);

	for (int i = 0; i < N; i++) {
		x = (i + 0.5) * dx;
		printf("%lf %lf %lf\n", x, uout[i], uexact[i]);
	}

	//hid_t        file_id, dataset_id, group_id, dataspace_id;  /* identifiers */
	//hsize_t      dims[1];
	//herr_t       status;
	//double* vec1 = (double*)malloc(N * sizeof(double));
	//free(vec1);

	///* HDF5: Create a new file to store velocity datasets. */
	//file_id = H5Fcreate("heat_trans_data.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//dims[0] = N;

	//dataspace_id = H5Screate_simple(1, dims, NULL);

	///* Create two groups to store initial values and output values in the file. */
	//group_id = H5Gcreate2(file_id, "/output", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	///* Create the datasets. */
	//dataset_id = H5Dcreate2(file_id, "/output/uout", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	////  printf("original dset_data[0]:%2d\n", dset_data[0]);

	// /* Write the first dataset. */
	//status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uout);

	///* Close the data space for the first dataset. */
	//status = H5Sclose(dataspace_id);

	///* Close the first dataset. */
	//status = H5Dclose(dataset_id);
	///* Close the group. */
	//status = H5Gclose(group_id);
	//status = H5Fclose(file_id);
	status = H5Sclose(dataspace_id1);
	status = H5Sclose(dataspace_id2);

	/* Close the first dataset. */
	status = H5Dclose(dataset_id1);
	status = H5Dclose(dataset_id2);
	/* Close the group. */
	status = H5Gclose(group_id);
	status = H5Fclose(file_id);

	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = VecDestroy(&u); CHKERRQ(ierr);
	ierr = VecDestroy(&uold); CHKERRQ(ierr);
	ierr = VecDestroy(&f); CHKERRQ(ierr);
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
}


// eulersch function
//double forward(double ujk, double ujmk, double ujpk, double ujkm, double ujkp, double x, double y) {
//	double u1, u2;
//	double ujknew;
//	double dt;
//	double f;
//
//	dt = Tend / (double)nsteps;
//	double dx;
//	dx = (double)L / (double)N;
//	double dy;
//	dy = (double)L / (double)N;
//	double pi = 4.0 * atan(1.0);
//
//	double CFL1 = kcond * dt / (dx * dx);
//	double CFL2 = kcond * dt / (dy * dy);
//	f = sin(L * pi * (x + y));
//	u1 = CFL1 * (ujpk - 2 * ujk + ujmk);
//	u2 = CFL2 * (ujkp - 2 * ujk + ujkm);
//	ujknew = ujk + 1 / (rho * c) * (f + u1 + u2);
//	return ujknew;
//}
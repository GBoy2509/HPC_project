// Finite-difference method to slove 2D heat transfer problem
// Square Omega of dimension 1 by 1
// Implicit Euler scheme for the time discretization
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
#define FILE2 "heat_trans_data.h5"

static char help[] = "Explicit EULER method\n\n";

int main(int argc, char** args) {
	// allocate two dimensional array (N+1)*(N+1)
	double** u = (double**)malloc(sizeof(float*) * (N + 2));
	for (int i = 0; i < (N + 2); i++) {
		u[i] = (double*)malloc(sizeof(double) * (N + 2));
	}
	double** unew = (double**)malloc(sizeof(float*) * (N + 2));
	for (int i = 0; i < (N + 2); i++) {
		unew[i] = (double*)malloc(sizeof(double) * (N + 2));
	}
	double** f = (double**)malloc(sizeof(float*) * (N));
	for (int i = 0; i < (N); i++) {
		f[i] = (double*)malloc(sizeof(double) * (N));
	}
	double x, y;
	double dx;
	dx = (double)L / (double)N;
	double dy;
	dy = (double)L / (double)N;
	FILE* init;
	FILE* uout;
	double pi = 4.0 * atan(1.0);
	printf("%d, %f, %f, %f\n", N, dx, dy, dt);

	double forward(double ujk, double ujmk, double ujpk, double ujkm, double ujkp, double x, double y);

	u[0][0] = 0.0;
	unew[0][0] = 0.0;
	u[N + 1][N + 1] = 0.0;
	unew[N + 1][N + 1] = 0.0;
	u[0][N + 1] = 0.0;
	unew[0][N + 1] = 0.0;
	u[N + 1][0] = 0.0;
	unew[N + 1][0] = 0.0;
	// grid init (x-d first, then y-d)
	// use initial condition
	for (int j = 1; j < (N + 1); j++) {
		y = (j - 0.5) * dy;
		for (int i = 1; i < (N + 1); i++) {
			x = (i - 0.5) * dx;
			u[i][j] = exp(x + y);
		}
	}
	// use boundary conditions
	// just set top & bottom as tau_g, left & right as tau_h
	for (int i = 1; i < (N + 1); i++) {
		u[i][0] = g;
		u[i][N + 1] = g;
		u[N + 1][i] = g;
		u[0][i] = g;

	}
	// for (int i = 1; i < (N + 1); i++) {
	//	 u[N + 1][i] = u[N][i] + dx * h1 / kcond;
	//	 u[0][i] = u[1][i] + dx * h2 / kcond;
	// }

	// test code
	init = fopen("init.csv", "w+");
	for (int j = 0; j < (N + 2); j++) {
		y = (j - 0.5) * dy;
		for (int i = 0; i < (N + 2); i++) {
			x = (i - 0.5) * dx;
			fprintf(init, "%f,%f,%f\n", x, y, u[i][j]);
		}
	}
	fclose(init);


	// Petsc Part
	MPI_Comm       comm;
	PetscMPIInt    rank;

	Mat            A, U, Uold, F;/* linear system matrix */
	Vec			   Ui, Uoldi;

	PetscReal      err, fvalue, value_to_set, norm_uold = 0.0, norm_u = 1.0;  /* norm of solution error */
	PetscErrorCode ierr;
	PetscInt       cc, col[3], num[N], maxit = 5000, rstart, rend, rlocal, clocal;
	PetscInt       max = 50000;
	PetscInt       n = N;
	PetscScalar    a = 1.0, value[3], vara[N];
	KSP		ksp;

	ierr = PetscInitialize(&argc, &args, (char*)0, help); if (ierr) return ierr;
	for (int i = 0; i < N; i++) {
		num[i] = i;
	}
	// MPI initial
	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);

	/*
	Create matrix U.  When using MatCreate(), the matrix format can
	be specified at runtime.
	Performance tuning note:  For problems of substantial size,
	preallocation of matrix memory is crucial for attaining good
	performance. See the matrix chapter of the users manual for details.
	We pass in nlocal as the "local" size of the matrix to force it
	to have the same parallel layout as the vector created above.
	*/
	ierr = MatCreate(PETSC_COMM_WORLD, &U); CHKERRQ(ierr);
	ierr = MatSetSizes(U, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
	ierr = MatSetFromOptions(U); CHKERRQ(ierr);
	ierr = MatSetUp(U); CHKERRQ(ierr);
	/* Identify the starting and ending mesh points on each
	processor for the interior part of the mesh. We let PETSc decide
	above. */
	ierr = MatGetOwnershipRange(U, &rstart, &rend); CHKERRQ(ierr);
	ierr = MatGetLocalSize(U, &rlocal, &clocal); CHKERRQ(ierr);

	for (int j = 1; j < (N + 1); j++) {
		for (int i = 1; i < (N + 1); i++) {
			ierr = MatSetValue(U, i, j, u[i][j], INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	/* Assemble the matrix */
	ierr = MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatView(U, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

	// Create coefs matrix A.
	ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
	ierr = MatSetSizes(A, rlocal, clocal, n, n); CHKERRQ(ierr);
	ierr = MatSetFromOptions(A); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(A, 3, NULL, 2, NULL);
	ierr = MatSetUp(A); CHKERRQ(ierr);
	// Set values to A
	MatGetOwnershipRange(A, &rstart, &rend);
	double CFL = dt * kcond / dx / dx / rho / c + dt * kcond / dy / dy / rho / c;
	int ii;
	if (!rstart)
	{
		rstart = 1;
		ii = 0; col[0] = 0; col[1] = 1; value[0] = 1.0 + 2.0 * CFL; value[1] = (- 1) * CFL;
		ierr = MatSetValues(A, 1, &ii, 2, col, value, INSERT_VALUES); CHKERRQ(ierr);
	}
	if (rend == n)
	{
		rend = n - 1;
		ii = n - 1; col[0] = n - 2; col[1] = n - 1; value[0] = (-1) * CFL; value[1] = 1.0 + 2.0 * CFL;
		ierr = MatSetValues(A, 1, &ii, 2, col, value, INSERT_VALUES); CHKERRQ(ierr);
	}
	value[0] = (-1) * CFL;
	value[1] = 1.0 + 2.0 * CFL;
	value[2] = (-1) * CFL;
	for (int i = rstart; i < rend; i++) {
		col[0] = i - 1; col[1] = i; col[2] = i + 1;
		MatSetValues(A, 1, &i, 3, col, value, INSERT_VALUES);
	}
	/* Assemble the matrix */
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	MatView(A, PETSC_VIEWER_STDOUT_WORLD);

	/* Duplicate the pattern of u to other vectors */
	MatDuplicate(U, &Uold);
	MatDuplicate(U, &F);

	/* Set value for heat supply, F */
	for (int j = 1; j < (N + 1); j++) {
		y = (j - 0.5) * dy;
		for (int i = 1; i < (N + 1); i++) {
			x = (i - 0.5) * dx;
			f[i][j] = sin(L * pi * (x + y)) * dt / rho / c;
		}
	}
	MatGetOwnershipRange(F, &rstart, &rend);
	for (int j = 1; j < (N + 1); j++) {
		for (int i = 1; i < (N + 1); i++) {
			ierr = MatSetValue(F, i, j, f[i][j], INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	MatAssemblyBegin(F);
	MatAssemblyEnd(F);
	MatView(F, PETSC_VIEWER_STDOUT_WORLD);

	VecCreate(PETSC_COMM_WORLD, &Ui);
	VecSetSizes(Ui, PETSC_DECIDE, n);
	VecSetFromOptions(Ui);
	VecGetOwnershipRange(Ui, &rstart, &rend);
	for (int i = rstart; i < rend; i++) {
		fvalue = 0.0;
		VecSetValue(Ui, i, fvalue, INSERT_VALUES);
	}
	VecAssemblyBegin(Ui);
	VecAssemblyEnd(Ui);
	VecView(Ui, PETSC_VIEWER_STDOUT_WORLD);
	VecDuplicate(Ui, &Uoldi);

	// Create the linear solver
	KSPCreate(PETSC_COMM_WORLD, &ksp);

	// calculate next step
	// backward euler scheme
	for (int iter = 0; iter < maxit; iter++) {
		MatCopy(U, Uold);
		MatAXPY(Uold, a, F, SAME_NONZERO_PATTERN);

		for (int i = 1; i < (n + 1); i++) {
			MatGetValues(Uold, n, num, 1, &i, vara);
			VecGetOwnershipRange(Uoldi, &rstart, &rend);
			for (int j = rstart; j < rend; j++) {
				value_to_set = vara[j];
				VecSetValue(Uoldi, j, value_to_set, INSERT_VALUES);
			}
			VecAssemblyBegin(Uoldi);
			VecAssemblyEnd(Uoldi);

			KSPSetOperators(ksp, A, A);
			KSPSetTolerances(ksp, 1.e-2 / n, 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);
			KSPSetFromOptions(ksp);
			KSPSolve(ksp, Uoldi, Ui);

			VecGetValues(Ui, n, num, vara);
			for (int j = 0; j < n; j++) {
				value_to_set = vara[j];
				MatSetValue(U, j, i, value_to_set, INSERT_VALUES);
			}
			MatAssemblyBegin(Ui, MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(Ui, MAT_FINAL_ASSEMBLY);
		}
		MatNorm(U, NORM_1, &norm_u);
		err = fabs(norm_u - norm_uold);
		if (err < 1.e-8) {
			PetscPrintf(comm, "End at %Dth iter, err=%g\n", iter + 1, err);
			break;
		}
		norm_uold = norm_u;
	}
	MatView(U, PETSC_VIEWER_STDOUT_WORLD);


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
	for (int i = 0; i < (N + 2); i++) {
		free(u[i]);
		free(unew[i]);
	}
	for (int i = 0; i < (N); i++) {
		free(f[i]);
	}
	free(u);
	free(unew);
	free(f);
	MatDestroy(&A);
	MatDestroy(&U);
	MatDestroy(&F);
	MatDestroy(&Uold);
	PetscFinalize();

	/* HDF5 initialization */
	hid_t        file_id, dataset_id, group_id, dataspace_id;  /* identifiers */
	hsize_t      dims[2];
	herr_t       status;
	int* vec1 = (int*)malloc((N + 2) * (N + 2) * sizeof(int));
	free(vec1);

	/* HDF5: Create a new file to store velocity datasets. */
	file_id = H5Fcreate("heat_trans_data.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = N + 2;
	dims[1] = N + 2;

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
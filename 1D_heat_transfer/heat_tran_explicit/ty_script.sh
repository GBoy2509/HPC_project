#!/bin/bash
#BSUB -J petsc-test
#BSUB -q ser
#BSUB -n 1
#BSUB -e %J-petsc.err
#BSUB -o %J-petsc.out

module purge
module load intel/2018.4
module load mpi/intel/2018.4

mpirun -np 1 ./heat_transfer_implicit.out \
  -log_view > $LSB_JOBID.log 2>&1

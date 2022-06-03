#!/bin/bash

#BSUB -J bash
#BSUB -q ser
#BSUB -n 1
#BSUB -e %bash.err
#BSUB -o %bash.out
#BSUB -R "span[ptile=1]"

module load intel/2018.4
module load gcc/9.3.0

./h5_combine


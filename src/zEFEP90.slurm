#!/bin/bash

#SBATCH --job-name=EFEP90
#SBATCH --partition=has28
#SBATCH -n 16
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=wangshenpeng@sjtu.edu.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=00:10:00

module load intel
make clean
make 
mpirun -n 16 ./a.exe

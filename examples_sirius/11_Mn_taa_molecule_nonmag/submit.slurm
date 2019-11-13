#!/bin/bash -l
#SBATCH --constraint=haswell
#SBATCH --partition=debug
#SBATCH --account=m3346
#SBATCH --time=00:30:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=1stbloodmpi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=longzhang2008@gmail.com

SECONDS=0                    # timing  

##export OMP_NUM_THREADS=4
##export MKL_NUM_THREADS=4
##srun -n 36 -c 4 --hint=nomultithread --unbuffered ./sirius.scf --mpi_grid="6 6"
##srun -n 36 -c 4 --hint=nomultithread --unbuffered ./sirius.scf --mpi_grid="2 2" --task=2 

export OMP_NUM_THREADS=16
ulimit -c unlimited
srun -n 16 ./elk --mpi_grid="4 4" > log_scf.txt
##srun -n 256 ./sirius.scf --task=2 --mpi_grid="16 16" > log_bands.txt

duration=$SECONDS            # timing

echo " "
echo " ================================================= "
echo " The job is done in time.                          "
echo " Time used for execution: $(($duration / 60)) min $(($duration % 60)) sec "
echo " ================================================= "
echo " " 

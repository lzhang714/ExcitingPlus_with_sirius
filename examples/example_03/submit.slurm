#!/bin/bash -l
#SBATCH --job-name=nioEPS
#SBATCH --constraint=haswell
#SBATCH --time=00:10:00 
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --partition=debug
#SBATCH -A mp261
#SBATCH --mail-type=ALL
#SBATCH --mail-user=longzhang2008@gmail.com

EXE=./elk 

SECONDS=0                    # for counting time  

srun -n 1 $EXE > elk.out
 
duration=$SECONDS            # for counting time

echo " "
echo " ================================================= "
echo " The job is done in time.                          "
echo " Time used for execution: $(($duration / 60)) min $(($duration % 60)) sec "
echo " ================================================= "
echo " " 

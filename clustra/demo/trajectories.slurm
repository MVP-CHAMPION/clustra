#!/bin/bash
#SBATCH -J trajectories
#SBATCH -A ccsd-burst
#SBATCH -D /home/ost
#SBATCH -t 01:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=4G
#SBATCH -e trajectories.e
#SBATCH -o trajectories.o
#SBATCH --mail-user=ost@ornl.gov
#SBATCH --mail-type=FAIL

hostname
cd ~/R/trajectories
pwd

## modules can vary on platforms
module load PE-gnu/3.0
module load R/3.6.0
echo "loaded R"
module list

## prevent warning when fork is used with MPI
export OMPI_MCA_mpi_warn_on_fork=0

Rscript R/trajectories.R

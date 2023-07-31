#!/bin/bash

#SBATCH -t time
#SBATCH -N nodes
#SBATCH -n cores
#SBATCH -p queue
#SBATCH -A allocation
#SBATCH --mail-user=email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=MOLECULE_NAME_CONFORMER_NUMBER

module reset
module load gaussian/09.e01

cd $SLURM_SUBMIT_DIR
pwd

#source $GAUSSIAN_DIR/bsd/g09.profile
export OMP_NUM_THREADS=$SLURM_NTASKS
export GAUSS_SCRDIR=$TMPDIR

echo "==================================="
echo "Running G09 on node:"
hostname
echo "===================================" 

g09 < input.dat > output.log

echo "=================================="
echo "Exiting"
echo "=================================="

exit;


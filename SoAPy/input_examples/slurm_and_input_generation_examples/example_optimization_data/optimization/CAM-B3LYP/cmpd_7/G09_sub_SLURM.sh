#!/bin/bash

#SBATCH -t 1-00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal_q
#SBATCH -A your_allocation
#SBATCH --mail-user=your_email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=H2O_cmpd_7_optimization_spherical_shell_CAM-B3LYP_STO-3G_1.0_10

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


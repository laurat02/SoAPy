#!/bin/bash

#SBATCH -t 6-00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p normal_q
#SBATCH -A crawdad
#SBATCH --mail-user=bshumberger@vt.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=B_glucose_cmpd_6_VCD_molecular_shell_CAM-B3LYP_cc-pVDZ_6-31G(d)_2.5_10

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


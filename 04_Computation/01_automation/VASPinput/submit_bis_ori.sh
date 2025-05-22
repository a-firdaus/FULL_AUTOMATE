#!/bin/sh

#SBATCH --job-name=
#SBATCH --time=hh:mm:ss
#SBATCH --mail-type=
#SBATCH --mail-user=
#SBATCH --ntasks=128
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32


module purge
module load VASP/5.4.4-intel-2019b # the one I've been using since I started in 2023
#module load VASP/6.3.2-gomkl-2022a
module list


cd $SLURM_SUBMIT_DIR

free -m
date
touch CHGCAR WAVECAR


srun vasp_std > $SLURM_SUBMIT_DIR/tempout
#mpirun sometimes works better


free -m
date
exit 0

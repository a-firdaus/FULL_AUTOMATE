#!/bin/sh

#SBATCH --job-name=aimd_nmc
#SBATCH --time=72:00:00
#SBATCH --mail-type=all 
#SBATCH --mail-user=azka.savanti.firdaus@vub.be
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

--output out.log
--error err.log

free -m
date
exit 0

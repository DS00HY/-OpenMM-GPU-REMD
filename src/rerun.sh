#!/bin/sh
#SBATCH -J re
#SBATCH -o re-%j.log
#SBATCH -e re-%j.err
#SBATCH -p baode_gpu_short
#SBATCH -t 5:00
#SBATCH -N 1
#SBATCH --gres=gpu:1
###SBATCH --ntasks-per-node=4
###SBATCH -N 1 --ntasks-per-node=6

echo begin re 1
mpirun -n 4 python restartREMD.py -debug
echo ok
####srun -n $SLURM_NTASKS python mysim.py >runrettwocg.log 2>&1
##mpirun -n 6 python mysim.py  >runrettwocg.log 2>&1
###bash run.sh
###python hh.py

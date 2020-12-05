#!/bin/sh
#SBATCH -J job_name
#SBATCH -o test-%j.log
#SBATCH -e test-%j.err
#SBATCH -p Sugon_40_long
#SBATCH -t 300
#SBATCH -N 1
###SBATCH --ntasks-per-node=4
###SBATCH -N 1 --ntasks-per-node=6


mpirun -n 4 python test.py
####srun -n $SLURM_NTASKS python mysim.py >runrettwocg.log 2>&1
##mpirun -n 6 python mysim.py  >runrettwocg.log 2>&1
###bash run.sh
###python hh.py

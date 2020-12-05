#!/bin/sh
#SBATCH -J job_name
#SBATCH -o job-%j.log
#SBATCH -e job-%j.err
#SBATCH -p baode_gpu_short
#SBATCH -t 200
#SBATCH -N 1
#SBATCH --gres=gpu:1
###SBATCH --ntasks-per-node=4
###SBATCH -N 1 --ntasks-per-node=6

echo begin
mpirun -n 4 python testREMD.py -debug
echo ok
####srun -n $SLURM_NTASKS python mysim.py >runrettwocg.log 2>&1
##mpirun -n 6 python mysim.py  >runrettwocg.log 2>&1
###bash run.sh
###python hh.py

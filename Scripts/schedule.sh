#!/bin/bash
##### Amount of cores per task
#SBATCH --cpus-per-task=32
##### Partition name
#SBATCH -p cpu
##### Name of job in queuing system
#SBATCH --job-name=QD
##### name of output file
#SBATCH --output="output.out"
##### name of error file
#SBATCH --error="error.err"

srun /home/czarnecki/LaoStoQuantumDots/bin/lao_sto_qd.x
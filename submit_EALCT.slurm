#!/bin/sh
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8192
#SBATCH --job-name=Fig_8_SB_n_AE
#SBATCH --partition=gpu --gres=gpu --constraint='gpu_32gb&gpu_v100'
#SBATCH --error=./job.%J.err
#SBATCH --output=./job.%J.out

module load matlab/r2021b
matlab -nodisplay -r "Fig_8_SB_n_AE, quit"
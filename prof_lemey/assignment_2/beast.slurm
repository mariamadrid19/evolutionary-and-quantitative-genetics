#!/bin/bash -l
#SBATCH --clusters=genius
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=00:05:00
#SBATCH --account=lp_edu_evol_quant_genetics 

cd $SLURM_SUBMIT_DIR
module purge
module use /apps/leuven/skylake/2021a/modules/all
module load beagle-lib
module load Beast/1.10.4-intel-2018a

beast YFV.xml

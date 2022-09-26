#! /bin/bash


######## login 
#SBATCH --job-name=main18_test
#SBATCH --output=./job-outs/main18_3.out
#SBATCH --error=./job-outs/main18_3.err


#SBATCH --account=pi-lhansen
#SBATCH --partition=standard
#SBATCH --cpus-per-task=3
#SBATCH --mem=36G
#SBATCH --time=7-00:00:00

####### load modules
module load python/booth/3.8/3.8.5  gcc/9.2.0

echo "$SLURM_JOB_NAME"

echo "Program starts $(date)"

python3 /home/bcheng4/Aiyagari/EGM_2Asset/main18.py 

echo "Program ends $(date)"


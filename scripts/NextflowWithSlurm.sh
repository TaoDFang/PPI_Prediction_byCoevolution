#!/usr/bin/env bash
#SBATCH --job-name "NextflowWithSlurm"
#SBATCH --output /shares/von-mering.imls.uzh/tao/nextflow/slurm_reports/output_%A_%a.txt # need to created own slurm_reports
#SBATCH --error /shares/von-mering.imls.uzh/tao/nextflow/slurm_reports/error_%A_%a.txt
#SBATCH --time=24:00:00 #96:00:00    #48:00:00
#SBATCH --mem=30G   # 5G here resource is to run nextlfow, not for acturally nextflow process ?
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=tao.fang@uzh.ch
#SBATCH --mail-type=ALL



# Use a conda environment where you have installed Nextflow
# (may not be needed if you have installed it in a different way)
# conda activate /data/tfang/conda-envs/nf-training
module load anaconda3
source activate /data/tfang/conda-envs/nf-training
cd /home/tfang/PPI_Prediction_byCoevolution/scripts
pwd
echo $CONDA_DEFAULT_ENV
nextflow run Query_coevolutionComputation_workflow.nf -c nextflow.config -profile slurm  -resume



# then run it by 
# chmod +x /home/tfang/PPI_Prediction_byCoevolution/scripts/NextflowWithSlurm.sh
# sbatch  /home/tfang/PPI_Prediction_byCoevolution/scripts/NextflowWithSlurm.sh
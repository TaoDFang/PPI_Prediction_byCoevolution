# PPI_Prediction_byCoevolution


## enviroment setup 
To run this repository, Nextflow and Singularity need to be insalled.

Nextflow installtion, follow official documentation:
https://www.nextflow.io/docs/latest/getstarted.html
or simpley done via conda install:  conda install -c bioconda nextflow=23.04.1
https://anaconda.org/bioconda/nextflow


Singularity installtaion follow official documentation:
https://docs.sylabs.io/guides/latest/user-guide/quick_start.html
or simply via conda install (Havent test this myself): conda install -c conda-forge singularity=3.8.7
https://anaconda.org/conda-forge/singularity

## Raw data and result computation 
To download all the raw data and to  generate all the data results needed for this project.

got to folder  PPI_Prediction_byCoevolution/scripts
then run one of the following:
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf  -c nextflow.config -profile singularity   -resume (on local machine)
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf  -c nextflow.config -profile slurm_withSingularity  -resume (On HPC with slurm)

## Paper figures
Download singularity container for this : 
singularity pull --arch amd64 library://tfang/base/py37_notebook:latest

The data needed to run notebooks can either be from last step that are generated from from scratch (This could take monthes depeonding on the avaiable computational resource)
Or The final cached results can be download from Zenodo at: 


to run notebooks:
go to folder PPI_Prediction_byCoevolution/notebooks and start the singularity container (by: singularity shell )
then inside container run:  jupyter notebook --no-browser --port=8036
check how to change or remove password here password here 


# PPI_Prediction_byCoevolution

in the end, copy paste to remmove all spelling error 

## enviroment setup 
To run this repository, Nextflow and Singularity need to be insalled.

Nextflow installtion:
[Official documentation](https://www.nextflow.io/docs/latest/getstarted.html) \
or simpley done via [conda install](https://anaconda.org/bioconda/nextflow):  
```
conda install -c bioconda nextflow=23.04.1
conda create -n py_nextflow --channel bioconda python=3.8 nextflow=23.04.1
```

The simplese strategy to use all necessery softwares and libraries for this project is to use the singularity container that we have built
Singularity installtaion:
[Official documentation](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html) \
or simply via [conda install](https://anaconda.org/conda-forge/singularity):  
```
conda install -c conda-forge singularity=3.8.7
conda create -n py_singularity --channel conda-forge python=3.8 singularity=3.8.7
```

the problem with conda installation, is that when then singulartyi is avaiable
mention here singularity is not necesseary 
(Havent test this myself)

If one has probelm to install singularity, all necessery softwares and library can be installed via conda. \
Check more details in documentation "containers/conda_envs/conda_installation.md" in this repository .


## Raw data and result computation 
mention time and size 

To download all the raw data and to  generate all the data results needed for this project.

got to folder  PPI_Prediction_byCoevolution/scripts \
then run one of the following: 
```
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf  -c nextflow.config -profile singularity   -resume (on local machine) \
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf  -c nextflow.config -profile slurm_withSingularity  -resume (On HPC with slurm)
```
for the customised parameters,you could directly modify thier values in configuration file "scripts/nextflow.config" \
or via nextflow the command line (e.g. --root_folder= "path to the location where you want to save all data")

here add more paramerter setting from here and introduce some pamameters , or direct add as comment line "scripts/nextflow.config" 


mention colabfold here

## Paper figures
Download singularity container for this :  
```
singularity pull --arch amd64 library://tfang/base/py37_notebook:latest
```

The data needed to run notebooks can either be from last step that are generated from from scratch (This could take monthes depeonding on the avaiable computational resource) \
Or The final cached results can be download from Zenodo at: \

to run notebooks: \
go to folder PPI_Prediction_byCoevolution/notebooks and start the singularity container by: 
```
singularity shell py37_notebook.sif
```
then inside container run: 
```
jupyter notebook --no-browser --port=8036 
```
then the notebooks are accessiable at: http://localhost:8036/ \
Remember to set variable "notebookData_folder" to the path where you save the data 
check how to change or remove password here password here  



# Testing

nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/mnt/mnemo6/tao" --conda_envs_path "/mnt/mnemo5/tao/anaconda3/envs" -c nextflow.config -profile standard  -resume

nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/mnt/mnemo6/tao"  -c nextflow.config -profile singularity  -resume
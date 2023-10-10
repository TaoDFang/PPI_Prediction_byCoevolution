# PPI_Prediction_byCoevolution

in the end, copy paste to remmove all spelling error 

## enviroment setup 
To run this repository, Nextflow and Singularity need to be insalled.

Nextflow installtion:
[Official documentation](https://www.nextflow.io/docs/latest/getstarted.html) \
or simpley done via [conda install](https://anaconda.org/bioconda/nextflow):  
```
conda install -c bioconda nextflow=23.04.1
```

Singularity installtaion:
[Official documentation](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html) \
or simply via [conda install](https://anaconda.org/conda-forge/singularity):  
```
conda install -c conda-forge singularity=3.8.7
```
mention here singularity is not necesseary 
(Havent test this myself)

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


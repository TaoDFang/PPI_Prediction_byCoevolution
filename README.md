# PPI_Prediction_byCoevolution


## environment setup 
To run this repository, Nextflow and Singularity need to be installed.

Nextflow installation:
[Official documentation](https://www.nextflow.io/docs/latest/getstarted.html) \
or simply done via [conda install](https://anaconda.org/bioconda/nextflow):  
```
conda create -n py_nextflow --channel bioconda python=3.8 nextflow=23.04.1
```

The simplest strategy to use all necessary software and libraries for this project is to use the singularity container that we have built
Singularity installation:
[Official documentation](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html) \
or simply via [conda install](https://anaconda.org/conda-forge/singularity)(Currently it's only avaiable for Linux):  
```
conda activate py_nextflow 
conda install -c conda-forge singularity=3.8.6 
```


If one has problems installing Singularity, all necessary software and libraries can be installed via conda (on Linux). \
Check more details in the documentation "containers/conda_envs/conda_installation.md" in this repository.


## Raw data and result computation 
To download all the raw data, generate paired alignment data, and compute DCA results for all protein pairs.

go to folder  PPI_Prediction_byCoevolution/scripts \
then run one of the following: 
```
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/home/tao/Data"  -c nextflow.config -profile singularity   -resume (on local machine) \
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/home/tao/Data"  -c nextflow.config -profile slurm_withSingularity  -resume (On HPC with slurm)
```
for the customized parameters, you could directly modify their values in the configuration file "scripts/nextflow.config" \
or via the command line (e.g. --root_folder= "path to the location where you want to save all data")

Warning: The whole computation could take months depending on the available computational resources and all the final results take up around 16TB of disk space

To run Alphafold-Multimer for the selected protein pairs, we use our customized paired alignment data as input to  [ColabFold (v1.3.0)](https://github.com/sokrypton/ColabFold/releases/tag/v1.3.0)

## Paper figures
Download the singularity container for this :  
```
singularity pull --arch amd64 library://tfang/base/py37_notebook:latest
```

The data needed to run notebooks can either be from the last step that is generated from scratch (This could take months depending on the available computational resource) \
Or The final cached results can be downloaded from Zenodo at: \

to run notebooks: \
go to folder PPI_Prediction_byCoevolution/notebooks and start the singularity container by: 
```
singularity shell py37_notebook.sif
```
then inside container run: 
```
jupyter notebook --no-browser --port=8036 
```
then the notebooks are accessible at http://localhost:8036/ \
Remember to set the variable "notebookData_folder" in the notebooks to the path where you save the data 




# Testing

test on VM (password is 123456)

nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/mnt/mnemo6/tao" --conda_envs_path "/mnt/mnemo5/tao/anaconda3/envs" -c nextflow.config -profile standard  -resume

nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/mnt/mnemo6/tao"  -c nextflow.config -profile singularity  -resume


nextflow run RawFastaFilesAndMetaData_workflow.nf -entry RawFastaFilesAndMetaData_workflow --root_folder "/Users/taof/Documents/PhD_Data"  -c nextflow.config -profile singularity  -resume


nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/Users/taof/Documents/PhD_Data" -c nextflow.config -profile singularity   -resume

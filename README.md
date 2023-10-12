# PPI_Prediction_byCoevolution

## environment setup 
To run this repository, Nextflow and Singularity need to be installed.

Nextflow installation:
[Official documentation](https://www.nextflow.io/docs/latest/getstarted.html) \
or simply done via [conda install](https://anaconda.org/bioconda/nextflow):  
```
conda create -n py_nextflow --channel bioconda python=3.8 nextflow=23.04.1
```

Singularity installation:
[Official documentation](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html) \
or simply via [conda install](https://anaconda.org/conda-forge/singularity)(Currently it's only avaiable for Linux):  
```
conda activate py_nextflow 
conda install -c conda-forge singularity=3.8.6 
```

The simplest strategy to use all necessary software and libraries for this project is to use the singularity containers that we have built. \
But if one has problems installing Singularity, all necessary software and libraries can be installed via conda (on Linux). \
Check more details in the documentation "containers/conda_envs/conda_installation.md" in this repository.


## Raw data and result computation 
To download all the raw data, generate paired alignment data, and compute DCA results for all protein pairs.

go to folder  PPI_Prediction_byCoevolution/scripts \
then run one of the following: 
```
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/home/tao"  -c nextflow.config -profile singularity   -resume (on the local machine) \
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/home/tao"  -c nextflow.config -profile slurm_withSingularity  -resume (On HPC with slurm)
nextflow run query2subject_homologousPPDetectionAndCompuation_workflow.nf --root_folder "/home/tao"  --conda_envs_path "/home/tao/anaconda3/envs" -c nextflow.config -profile standard   -resume (on the local machine when the singularity is not available) \
```
for the customized parameters, you could directly modify their values in the configuration file "scripts/nextflow.config" \
or via the command line (e.g. --root_folder= "path to the location where you want to save all data")

Warning: The whole computation could take months depending on the available computational resources and all the final results take up around 16TB of disk space \
This dataset is too large to share online so is only available upon request.


To run Alphafold-Multimer for the selected protein pairs, we use the generated customized paired alignment data as input to  [ColabFold (v1.3.0)](https://github.com/sokrypton/ColabFold/releases/tag/v1.3.0)

## Paper figures
Download the singularity container for this :  
```
singularity pull --arch amd64 library://tfang/base/py38_notebook:latest
```

The data needed to run notebooks can either be from the last step that is generated from scratch (This could take months depending on the available computational resource) \
Or The final cached results can be downloaded from Zenodo at: https://zenodo.org/record/8429824

to run notebooks: \
go to folder PPI_Prediction_byCoevolution/notebooks and start the singularity container by: 
```
singularity shell py38_notebook.sif
```
then inside container run: 
```
jupyter notebook --no-browser --port=8036 
```
then the notebooks are accessible at http://localhost:8036/ \
Remember to set the variable "notebookData_folder" in the notebooks to the location where you save the data 

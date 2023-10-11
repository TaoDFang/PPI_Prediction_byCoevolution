To install conda enviroments "py37_notebook" and "sequence_tools_environment", simply run following code : \
```
conda env  create  -f ~/PPI_Prediction_byCoevolution/containers/conda_envs/sequence_tools_environment.yml  \
conda env  create  -f ~/PPI_Prediction_byCoevolution/containers/conda_envs/py37_notebook_environment.yml 
```



To install conda environment "py37_pydca" with python library [pydca](https://github.com/KIT-MBS/pydca) installed. \
As the NumPy package always causes conflicts for this library, we recommend installing it through the following commands : 
```

conda create -n py37_pydca python==3.7
conda activate py37_pydca

pip install pydca

pip uninstall -y numpy 
pip uninstall -y numba

conda install -c numba numba=0.46.0

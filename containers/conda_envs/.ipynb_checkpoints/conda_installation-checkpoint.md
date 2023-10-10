To install conda enviroments "py37_notebook" and "sequence_tools_environment", simply run following code : \
```
conda env  create  -f ~/PPI_Prediction_byCoevolution/containers/conda_envs/sequence_tools_environment.yml  \
conda env  create  -f ~/PPI_Prediction_byCoevolution/containers/conda_envs/py37_notebook_environment.yml 
```



To install conda enviroment "py37_pydca" with python library [pydca](https://github.com/KIT-MBS/pydca) installed. \
As the numbpy package always cause the conflicts for this library,  so we recommend to install it through following commands : 
```

conda create -n py37_pydca python==3.7
conda activate py37_pydca

pip install pydca

pip uninstall -y numpy 
pip uninstall -y numba

conda install -c numba numba=0.46.0
```
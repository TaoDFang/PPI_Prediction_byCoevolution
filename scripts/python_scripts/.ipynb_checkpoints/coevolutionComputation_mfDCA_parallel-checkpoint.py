import argparse
import os 
# import pandas as pd , py37_pydca envs dont have pandas package to mess up the numpy version 
import csv
from datetime import datetime



from DCA_computation import pydca_mfdca_FN_compresse


#this is to force pydca use limit number of cpus per multiprocessing  process, its also works in nextflow pipeline, 
os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4

os.environ["MKL_NUM_THREADS"] = "4" # export MKL_NUM_THREADS=6

os.environ["NUMBA_NUM_THREADS"] = "4" # export NUMBA_NUM_THREADS=6



## to run pydca, its import to make sure these three package are in correct version 
## use py37_pydca enviroment 
# (py37_pydca) tao@deimos:~$ conda remove numpy
# (py37_pydca) tao@deimos:~$ pip uninstall numpy 
# (py37_pydca) tao@deimos:~$ conda install numpy=1.15.4
# np.__version__, numba.__version__, llvmlite.__version__

# (py37_pydca) tao@gaia:~$ pip  list | grep  numpy
# numpy                         1.15.4
# (py37_pydca) tao@gaia:~$ conda   list | grep  numpy
# numpy                     1.15.4          py37h8b7e671_1002    conda-forge

# version works 
# ('1.15.4', '0.46.0', '0.30.0')
# or ('1.18.4', '0.46.0', '0.30.0')
# or ('1.19.2', '0.46.0', '0.30.0')

# choosed version 
# ('1.18.4', '0.46.0', '0.30.0')
import llvmlite
import numba
import numpy as np 

print("np.__version__, numba.__version__, llvmlite.__version__:",np.__version__, numba.__version__, llvmlite.__version__)


if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Compute_allPPI_testOneIndex')
    parser.add_argument('-dipath','--IndexDCA_coevolutoin_path', type=str, help='IndexDCA_coevolutoin_path')
    parser.add_argument('-i','--idx', type=str, help='idx to files to be computed')
    

    args = parser.parse_args()
    IndexDCA_coevolutoin_path=args.IndexDCA_coevolutoin_path
    index_count=int(args.idx)




    # block_currentSpe_allPPIs_pps_forDCA_frame=pd.read_csv(IndexDCA_coevolutoin_path+str(index_count)+".csv",
    #                                                         header=None,index_col=None,sep="\t")

    # print("block_currentSpe_allPPIs_pps_forDCA_frame.shape",block_currentSpe_allPPIs_pps_forDCA_frame.shape)
    # block_currentSpe_allPPIs_pps_forCoevolution_list=block_currentSpe_allPPIs_pps_forDCA_frame.values.tolist()

    block_currentSpe_allPPIs_pps_forCoevolution_list=list()
    with open(IndexDCA_coevolutoin_path+str(index_count)+".csv", "r") as file:
        my_reader = csv.reader(file, delimiter="\t")
        # next(my_reader, None)  # skip the headers, no header in this file 
        for row in my_reader:
            block_currentSpe_allPPIs_pps_forCoevolution_list.append(row)
    print("len(block_currentSpe_allPPIs_pps_forCoevolution_list):",len(block_currentSpe_allPPIs_pps_forCoevolution_list))        


    now = datetime.now()
    print("Coevolution comp stat now =", now)

    runned_count=0
    for i in range(len(block_currentSpe_allPPIs_pps_forCoevolution_list)): 
        print(i)
        DCA_l=block_currentSpe_allPPIs_pps_forCoevolution_list[i]
        pydca_mfdca_FN_compresse(DCA_l)
        runned_count+=1

    now = datetime.now()
    print("Coevolution comp end now =", now)
    print("runned_count:",runned_count)









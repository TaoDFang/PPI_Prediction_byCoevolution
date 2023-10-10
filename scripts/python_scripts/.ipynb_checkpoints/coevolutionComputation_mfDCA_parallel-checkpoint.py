import argparse
import os 
import csv
from datetime import datetime



from DCA_computation import pydca_mfdca_FN_compresse


#this is to force pydca use limit number of cpus per multiprocessing  process, its also works in nextflow pipeline, 
os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4

os.environ["MKL_NUM_THREADS"] = "4" # export MKL_NUM_THREADS=6

os.environ["NUMBA_NUM_THREADS"] = "4" # export NUMBA_NUM_THREADS=6


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









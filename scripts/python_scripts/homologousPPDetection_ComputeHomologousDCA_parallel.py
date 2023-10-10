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

    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--Subject_tupleList', type=str, help='Subject_tupleList')
    parser.add_argument('-c','--CoEvo_data_folder', type=str, help='CoEvo_data_folder')
    parser.add_argument('-i','--idx', type=str, help='idx to files to be computed')
    
    
    args = parser.parse_args()
    Subject_tupleList=args.Subject_tupleList
    CoEvo_data_folder=args.CoEvo_data_folder
    index_count=int(args.idx)
    
    
    Subject_tupleList=Subject_tupleList.split("_")
    Subject_tupleList=[(Subject_tupleList[i],Subject_tupleList[i+1]) for i in range(0,len(Subject_tupleList)-1,2)]
    
    
    for currentSubject_EggNOG_maxLevel,currentSubject_TaxID in Subject_tupleList:
        print(currentSubject_EggNOG_maxLevel,currentSubject_TaxID,index_count)
                
        input_root_folder=CoEvo_data_folder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_eggNOGfilteredData/"
        IndexDCA_coevolutoin_path=input_root_folder+"coevolutoin_computation_IndexDCA/"
        print("IndexDCA_coevolutoin_path:",IndexDCA_coevolutoin_path)


        # block_currentSpe_allPPIs_pps_forDCA_frame=pd.read_csv(IndexDCA_coevolutoin_path+str(index_count)+".csv",
        #                                                         header=None,index_col=None,sep="\t")        
        # print("block_currentSpe_allPPIs_pps_forDCA_frame.shape",block_currentSpe_allPPIs_pps_forDCA_frame.shape)
        # block_currentSpe_allPPIs_pps_forCoevolution_list=block_currentSpe_allPPIs_pps_forDCA_frame.values.tolist()
        
        
        block_currentSpe_allPPIs_pps_forCoevolution_list=list()
        
        try:
            with open(IndexDCA_coevolutoin_path+str(index_count)+".csv", "r") as file:
                my_reader = csv.reader(file, delimiter="\t")
                # next(my_reader, None)  # skip the headers, no header in this file 
                for row in my_reader:
                    block_currentSpe_allPPIs_pps_forCoevolution_list.append(row)
            print("len(block_currentSpe_allPPIs_pps_forCoevolution_list):",len(block_currentSpe_allPPIs_pps_forCoevolution_list))  
        except:
            assert index_count>0
            print("It seems now only few samples wait to be compuated, index_count is :",str(index_count))
        

        now = datetime.now()
        print("Coevolution comp stat now =", now)

        runned_count=0
        for i in range(len(block_currentSpe_allPPIs_pps_forCoevolution_list)): 
        #for i in range(1): 
            print(i)
            DCA_l=block_currentSpe_allPPIs_pps_forCoevolution_list[i]
            pydca_mfdca_FN_compresse(DCA_l)
            runned_count+=1

        now = datetime.now()
        print("Coevolution comp end now =", now)
        print("runned_count:",runned_count)




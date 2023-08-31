import argparse
import os 
import glob
import pandas as pd 
import multiprocessing as mp
from multiprocessing import get_context
import pickle
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

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-dpath','--DCA_coevolutoin_path', type=str, help='DCA_coevolutoin_path')
    parser.add_argument('-m','--currentSpeMSAGapsFilteringMetaFolder', type=str, help='currentSpeMSAGapsFilteringMetaFolder')
    parser.add_argument('-acsv','--PPIInfoBeforeCoEvoComp_csv', type=str, help='PPIInfoBeforeCoEvoComp_csv')
    parser.add_argument('-nf90f','--pairedMSA_Nf90_folder', type=str, help='pairedMSA_Nf90_folder')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    
    
    args = parser.parse_args()
    currentSpeMSAGapsFilteringMetaFolder=args.currentSpeMSAGapsFilteringMetaFolder
    PPIInfoBeforeCoEvoComp_csv=args.PPIInfoBeforeCoEvoComp_csv
    DCA_coevolutoin_path=args.DCA_coevolutoin_path
    pairedMSA_Nf90_folder=args.pairedMSA_Nf90_folder
    mp_task_nums=int(args.mp_task_nums)

    
    with open(currentSpeMSAGapsFilteringMetaFolder+'fasta_protein_lens_dict.pickle', 'rb') as handle:
        fasta_protein_lens=pickle.load(handle)

    # use  postive and negative ppi to only use pp that have large Nf90 value, and after removing deep homologs  
    currentSpe_allPPIs_beforeCoEvoComp_frame = pd.read_csv(PPIInfoBeforeCoEvoComp_csv,header=0,index_col=None,sep="\t")
    print("currentSpe_allPPIs_beforeCoEvoComp_frame.shape:",currentSpe_allPPIs_beforeCoEvoComp_frame.shape)
    currentSpe_allPPIs_beforeCoEvoComp_info=currentSpe_allPPIs_beforeCoEvoComp_frame.values.tolist()
    currentSpe_allPPIs_pps=[(p1,p2) for p1, p2,_ in currentSpe_allPPIs_beforeCoEvoComp_info]
    print("len(currentSpe_allPPIs_pps):",len(currentSpe_allPPIs_pps))
    


    #if re-run, first check if file is already existed, need to change pythion code 
    existed_pydcaFNAPC_files=glob.glob(DCA_coevolutoin_path+"*_pydcaFNAPC_array.npz")
    print("len(existed_pydcaFNAPC_files):",len(existed_pydcaFNAPC_files))
    print("existed_pydcaFNAPC_files[0:3]:",existed_pydcaFNAPC_files[0:3])
    existed_pydcaFNAPC_files=[os.path.basename(f) for f in existed_pydcaFNAPC_files]
    existed_pydcaFNAPC_pps=[f.split("and") for f in existed_pydcaFNAPC_files]
    existed_pydcaFNAPC_pps=[(p1,p2[:-21]) for p1, p2 in existed_pydcaFNAPC_pps]
    existed_pydcaFNAPC_pp_dict=dict([(pp,1) for pp in existed_pydcaFNAPC_pps])



    #%%time 
    currentSpe_allPPIs_pps_forDCA=[pp for pp in currentSpe_allPPIs_pps if pp not in existed_pydcaFNAPC_pp_dict]
    currentSpe_allPPIs_pps_ArgForDCA=[(p1,p2,fasta_protein_lens[p1],fasta_protein_lens[p2],pairedMSA_Nf90_folder,DCA_coevolutoin_path) for p1 , p2 in currentSpe_allPPIs_pps_forDCA]
    #for dedub reson, delete it later 
    currentSpe_allPPIs_pps_ArgForDCA=currentSpe_allPPIs_pps_ArgForDCA[0:3]
    print("len(currentSpe_allPPIs_pps_ArgForDCA):",len(currentSpe_allPPIs_pps_ArgForDCA))


    if len(currentSpe_allPPIs_pps_ArgForDCA)>0:


        now = datetime.now()
        print("DCA comp stat now =", now)

        pool=get_context("spawn").Pool(mp_task_nums)#get_context("spawn").Pool(30) #mp.Pool(20)

        # Here we see  broken pipe problem 
        # more interesting thing is that when DCA compuation is still running , especial other pipeline for folowing 
        # speices alreays start. 
        # from loging file we can see MI part never stat . so the problem is cause by DCA part
        # so DCA part broken, current for loop for current speceis finished and pipelie for other speceis starts 
        # but in term of current for loop. part after DCA computation were not started

        dca_results=pool.map(pydca_mfdca_FN_compresse,currentSpe_allPPIs_pps_ArgForDCA)

        # time.sleep(10)# add this to prevent broken pipeline problme ??
        # # to prevent pool is close before  geting all return ???

        pool.close() 

    now = datetime.now()
    print("DCA comp end now =", now)




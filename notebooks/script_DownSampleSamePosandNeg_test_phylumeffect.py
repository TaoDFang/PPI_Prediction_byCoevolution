#!/usr/bin/env python
# coding: utf-8

# '''
# http://localhost:8206/lab/workspaces/auto-Y/tree/code/MNF/notebooks/Extend2OtherBacterialSpe/CoEvo_EggNOG_preprocessing.ipynb
# '''
# need enviroment py37_pydca

# In[ ]:



from __future__ import division
import argparse
import re
import multiprocessing as mp
from multiprocessing import get_context
import glob
import pandas as pd
import sys
import os
import math
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
from Bio import AlignIO

from shutil import copyfile
from collections import defaultdict
import numpy as np
import matplotlib as plt
import matplotlib.pyplot as plt
import csv
import gzip
import time
import logging
from datetime import datetime


import warnings
import pickle

import subprocess



import itertools

import smtplib
import ssl
import yaml


import random 
random.seed(10)



os.getcwd()





# In[ ]:


# set number of thread used by pydca 
# better set these before actuall load python library 
#https://stackoverflow.com/questions/30791550/limit-number-of-threads-in-numpy
    

import os
os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4

os.environ["MKL_NUM_THREADS"] = "4" # export MKL_NUM_THREADS=6

os.environ["NUMBA_NUM_THREADS"] = "4" # export NUMBA_NUM_THREADS=6
#set tocken


# In[ ]:


import sys
print(sys.executable)
print(sys.version)
print(sys.version_info)


# In[ ]:


from pydca.meanfield_dca import meanfield_dca
from pydca.plmdca import plmdca

import llvmlite
import numba
import numpy as np 


# In[ ]:


## ???!!! to run pydca, its import to make sure these three package are in correct version 
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

np.__version__, numba.__version__, llvmlite.__version__


# In[ ]:





# In[ ]:


#get_ipython().run_line_magic('reload_ext', 'autoreload')
#get_ipython().run_line_magic('autoreload', '2')

sys.path.append("/mnt/mnemo5/tao/code/MNF/src/tao_utilities")

from create_singleMSA import ParseCurSpeFastaByProteins_STRNG1105
from create_singleMSA import combinedRBH_STRING1105
from create_singleMSA import fun_newSingleMSA_EggNOG_OrthologousGroup_faidx
from create_singleMSA import fun_newSingleMSA_EggNOG_OrthologousGroup_phmmer
from create_singleMSA import ClustoMSA
from create_singleMSA import hmmbuild
from create_singleMSA import hmmalign
from create_singleMSA import removeGapsANDtrackAAPosOfSeqInMSA


from  create_pairedMSA import get_pairedMSA_inOneRun
from PPI_benchmark_filter import  get_PPIwithLimitedProFreByOr
from PPI_benchmark_filter import  getSameProteinRatio
from DCA_computation import pydca_mfdca_FN_compresse
from MI_computation import MI_calulation_new
from collect_topCoEvos import get_maxBetValue_dict_pydcaFNAPC_array_npz
from collect_topCoEvos import get_maxBetValue_dict_MI_apc_allResidues_npz
from collect_topCoEvos import getFrame2Dict_bet
from PPI_benchmark_combine import combineManyInfo_simpleVersion


# In[ ]:





# Here use ortholog group from eggNOG provided by Damian to  reconstruct profileHMM 
# The one thing should be notice here is that I should better use use orgholog group provided by Damian 
# instead of from eggNOG website because their one ortholog group could have multiple proteins(paralogs) 
# from  same specie http://eggnog5.embl.de/#/app/results.
# I think this is also the reason that eggNOG doesn't provide profileHMM for such orthologous groups . 
# Good news for me is that by only keeping  one paralog per-species, the number of protein sequence  
# in each orthologous group would be much reduced 
# 

# In[ ]:





# # Set some importtant path,folders and global variables 
# 

# In[ ]:





if __name__ == '__main__':


    # set email notification

    conf = yaml.safe_load(open('/mnt/mnemo5/tao/conf/application.yml'))
    port = 587  # For starttls
    smtp_server = "smtp.uzh.ch"
    sender_email = conf['user']['email']
    receiver_email = conf['user']['email']
    password =conf['user']['password']
    context = ssl.create_default_context()




    parser = argparse.ArgumentParser(description='script_DownsampleSamePosandNeg_test_phylumeffect')
    parser.add_argument('-l','--EggNOG_maxLevel', type=str, help='EggNOG_maxLevel')
    parser.add_argument('-i','--currentSpe_TaxID', type=str, help='currentSpe_TaxID')
    parser.add_argument('-s','--DownSample_strategy', type=str, help='DownSample_strategy')
    parser.add_argument('-m','--DownSample_size', type=int, help='DownSample_size')
    

    args = parser.parse_args()
    EggNOG_maxLevel=args.EggNOG_maxLevel
    currentSpe_TaxID=args.currentSpe_TaxID
    DownSample_strategy=args.DownSample_strategy
    DownSample_size=args.DownSample_size
    

    # # -----------------prepare benchmark dataset -----------------------------


    newSTRING_rootFolder="/mnt/mnemo6/tao/PPI_Coevolution/STRING_data_11.5/"    
    CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING11.5/"


    input_root_folder=CoEvo_data_folder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_eggNOGfilteredData/"
    pairedMSA_Nf90_folder = input_root_folder+"DownSample_"+DownSample_strategy+str(DownSample_size)+"_pair_MSA_Nf90_PasteAlign/"


    DCA_coevolutoin_path=input_root_folder+"DownSample_"+DownSample_strategy+str(DownSample_size)+"_coevolutoin_result_DCA/"
    MI_coevolutoin_path=input_root_folder+"DownSample_"+DownSample_strategy+str(DownSample_size)+"_coevolutoin_result_MI/"

    Benchmark_folder=input_root_folder+"DownSample_"+DownSample_strategy+str(DownSample_size)+"_testPhylaEffectDownSampleSamePPs_Benchmark/"
    
    print(pairedMSA_Nf90_folder,DCA_coevolutoin_path,MI_coevolutoin_path,Benchmark_folder)
    if not os.path.exists(DCA_coevolutoin_path):
        os.makedirs(DCA_coevolutoin_path)

    if not os.path.exists(MI_coevolutoin_path):
        os.makedirs(MI_coevolutoin_path)
        
    if not os.path.exists(Benchmark_folder):
        os.makedirs(Benchmark_folder)
        
    
    currentSpeMiddleDataPath=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_MiddleData/"
    with open(currentSpeMiddleDataPath+'fasta_protein_lens_dict.pickle', 'rb') as handle:
        fasta_protein_lens=pickle.load(handle)
                                       
    print(fasta_protein_lens["511145.b0014"],fasta_protein_lens["511145.b0055"],fasta_protein_lens["511145.b0095"],fasta_protein_lens["511145.b0436"])



    SamePosandNeg_Benchmark=CoEvo_data_folder+"AllLevel_testPhylaEffectSamePosandNeg_Benchmark/"

    SamePosandNeg_frame=pd.read_csv(SamePosandNeg_Benchmark+"SamePosandNeg_STRING115EggNog"+EggNOG_maxLevel+"Spe"+currentSpe_TaxID+".benchmark",
                                                                                                                header=0,index_col=None,sep="\t")
    print(SamePosandNeg_frame.shape)
    SamePosandNeg_list=SamePosandNeg_frame.loc[:,["STRING_ID1","STRING_ID2","len1","len2",]].values.tolist()
    print(len(SamePosandNeg_list))
    SamePosandNeg_frame.head(n=3)
    
    
    SamePosandNeg_PosPPs_list=SamePosandNeg_frame.loc[SamePosandNeg_frame["benchmark_status"]=="P",["STRING_ID1","STRING_ID2",]].values.tolist()
    SamePosandNeg_PosPPs_list=[tuple(pp) for pp in SamePosandNeg_PosPPs_list]
    SamePosandNeg_NegPPs_list=SamePosandNeg_frame.loc[SamePosandNeg_frame["benchmark_status"]=="N",["STRING_ID1","STRING_ID2",]].values.tolist()
    SamePosandNeg_NegPPs_list=[tuple(pp) for pp in SamePosandNeg_NegPPs_list]
    print(len(SamePosandNeg_PosPPs_list),len(SamePosandNeg_NegPPs_list),SamePosandNeg_PosPPs_list[0:3])
    
    #if re-run, first check if file is already existed, need to change pythion code 

    existed_pydcaFNAPC_files=glob.glob(DCA_coevolutoin_path+"*_pydcaFNAPC_array.npz")
    print("len(existed_pydcaFNAPC_files):",len(existed_pydcaFNAPC_files))
    print(existed_pydcaFNAPC_files[0:3])

    existed_pydcaFNAPC_files=[os.path.basename(f) for f in existed_pydcaFNAPC_files]
    existed_pydcaFNAPC_pps=[f.split("and") for f in existed_pydcaFNAPC_files]
    existed_pydcaFNAPC_pps=[(p1,p2[:-21]) for p1, p2 in existed_pydcaFNAPC_pps]

    existed_pydcaFNAPC_pp_dict=dict([(pp,1) for pp in existed_pydcaFNAPC_pps])
    
    #%%time 
    currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForDCA=[(p1,p2,l1,l2,pairedMSA_Nf90_folder,DCA_coevolutoin_path) for p1,p2,l1,l2 in SamePosandNeg_list if (p1,p2) not in existed_pydcaFNAPC_pp_dict]
    print("len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForDCA):",len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForDCA))
    #print(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForDCA[0])
    
    if len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForDCA)>0:


        now = datetime.now()
        print("DCA comp stat now =", now)

        pool=get_context("spawn").Pool(30)#get_context("spawn").Pool(30) #mp.Pool(20)

        # Here we see  broken pipe problem 
        # more interesting thing is that when DCA compuation is still running , especial other pipeline for folowing 
        # speices alreays start. 
        # from loging file we can see MI part never stat . so the problem is cause by DCA part
        # so DCA part broken, current for loop for current speceis finished and pipelie for other speceis starts 
        # but in term of current for loop. part after DCA computation were not started

        dca_results=pool.map(pydca_mfdca_FN_compresse,currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForDCA)

        time.sleep(10)# add this to prevent broken pipeline problme ??
        # to prevent pool is close before  geting all return ???

        pool.close() 




    now = datetime.now()
    print("DCA comp end now =", now)
    
    
    # ### MI

    # In[ ]:


    #if re-run, first check if file is already existed, need to change pythion code 

    existed_MIAPC_files=glob.glob(MI_coevolutoin_path+"*_apc_allResidues.npz")
    print("len(existed_MIAPC_files):",len(existed_MIAPC_files))
    print(existed_MIAPC_files[0:3])

    existed_MIAPC_files=[os.path.basename(f) for f in existed_MIAPC_files]
    existed_MIAPC_pps=[f.split("and") for f in existed_MIAPC_files]
    existed_MIAPC_pps=[(p1,p2[:-20]) for p1, p2 in existed_MIAPC_pps]

    existed_MIAPC_pp_dict=dict([(pp,1) for pp in existed_MIAPC_pps])
    
    
    #%%time 
    currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForMI=[(p1,p2,l1,l2,pairedMSA_Nf90_folder,MI_coevolutoin_path) for p1,p2,l1,l2 in SamePosandNeg_list if (p1,p2) not in existed_MIAPC_pp_dict]
    print("len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForMI):",len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForMI))
    #print(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForMI[0])
    
    

    if len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForMI)>0:
        now = datetime.now()
        print(" MI comp stat now =", now)

        #pool=get_context("spawn").Pool(50)#get_context("spawn").Pool(50) #mp.Pool(50)
        pool=get_context("spawn").Pool(50) # here set to 30 to see if we can solve broken pipe problem 
        # after 30 we still see problem 
        # more interesting thing is that when DCA compuation is still running ,  especial other pipeline for folowing 
        # speices alreays start. 
        # from loging file we can see MI part never stat . so the problem is cause by DCA part
        MI_results=pool.map(MI_calulation_new,currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForMI)

        time.sleep(10)# add this to prevent broken pipeline problme ??
        # to prevent pool is close before  geting all return ???

        pool.close() 


    now = datetime.now()
    print("MI comp end now =", now)
    
    currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForGetDCAMax=[(p1,p2,l1,l2,DCA_coevolutoin_path) for p1 , p2,l1,l2 in SamePosandNeg_list]
    print("len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForGetDCAMax):",len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForGetDCAMax))

    if not os.path.exists(Benchmark_folder+"max_pydcaFNAPC_frame.csv"):

        pool=mp.Pool(30) 
        max_pydcaFNAPC_list=pool.map(get_maxBetValue_dict_pydcaFNAPC_array_npz,currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForGetDCAMax)
        pool.close() 
        max_pydcaFNAPC_frame=pd.DataFrame(max_pydcaFNAPC_list,columns=["currentSpe_pro1", "currentSpe_pro2",
                                                                       "maxValue_all_idx_row","maxValue_all_idx_col","maxValue_all",
                                                                       "maxValue_bet_idx_row","maxValue_bet_idx_col","maxValue_bet"])
        max_pydcaFNAPC_frame.to_csv(Benchmark_folder+"max_pydcaFNAPC_frame.csv",header=True, index=None,sep="\t")


    currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForGetMIMax=[(p1,p2,l1,l2,MI_coevolutoin_path) for p1 , p2,l1,l2 in SamePosandNeg_list]
    print("len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForGetMIMax):",len(currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForGetMIMax))
    if not os.path.exists(Benchmark_folder+"max_MIApcAllResidues_frame.csv"):


        # why this section is slow then last section , because of MI has many more outpt in the fodler ???
        pool=mp.Pool(30) 
        max_MIApcAllResidues_list=pool.map(get_maxBetValue_dict_MI_apc_allResidues_npz,currentSpe_STRINGcurrenttestPhylaEffectSamePPI_allPPI_ArgForGetMIMax)
        pool.close() 
        max_MIApcAllResidues_frame=pd.DataFrame(max_MIApcAllResidues_list,columns=["currentSpe_pro1", "currentSpe_pro2",
                                                                       "maxValue_all_idx_row","maxValue_all_idx_col","maxValue_all",
                                                                       "maxValue_bet_idx_row","maxValue_bet_idx_col","maxValue_bet"])
        max_MIApcAllResidues_frame.to_csv(Benchmark_folder+"max_MIApcAllResidues_frame.csv",header=True, index=None,sep="\t")



    max_pydcaFNAPC_frame=pd.read_csv(Benchmark_folder+"max_pydcaFNAPC_frame.csv",
                                    header=0,index_col=None,sep="\t")
    print(max_pydcaFNAPC_frame.shape)
    max_pydcaFNAPC_bet_dict=getFrame2Dict_bet(max_pydcaFNAPC_frame)
    print(list(max_pydcaFNAPC_bet_dict.items())[0:3])


    # In[ ]:


    max_MIApcAllResidues_frame=pd.read_csv(Benchmark_folder+"max_MIApcAllResidues_frame.csv",
                                    header=0,index_col=None,sep="\t")
    print(max_MIApcAllResidues_frame.shape)
    max_MIApcAllResidues_bet_dict=getFrame2Dict_bet(max_MIApcAllResidues_frame)



    allPPI_allInfo_frame=combineManyInfo_simpleVersion(SamePosandNeg_PosPPs_list,SamePosandNeg_NegPPs_list,
                                                      fasta_protein_lens,
                                                      max_MIApcAllResidues_bet_dict,max_pydcaFNAPC_bet_dict,)
    allPPI_allInfo_frame.to_csv(Benchmark_folder+"allPPI_allInfo_frame.csv",header=True,index=None,sep="\t")


    # In[ ]:



    # In[ ]:


    allPPI_allInfo_frame=pd.read_csv(Benchmark_folder+"allPPI_allInfo_frame.csv",
                                     header=0,index_col=None,sep="\t")
    print(allPPI_allInfo_frame.shape)
    allPPI_allInfo_frame.head(n=3)



 
 





    print("nice, whole pipeline finished for  species and eggnog level : ",currentSpe_TaxID ,";",EggNOG_maxLevel)

    message = "nice, whole pipeline finished for  species and eggnog level : "+currentSpe_TaxID +";"+EggNOG_maxLevel
    with smtplib.SMTP(smtp_server, port) as server:
        server.ehlo()  # Can be omitted
        server.starttls(context=context)
        server.ehlo()  # Can be omitted
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_email, message)
    # except: 
    #     print("oops , some filter steps didnt pass, read log file; for  species and eggnog level:  ",currentSpe_TaxID ,";",EggNOG_maxLevel)
    #     message="oops , some filter steps didnt pass, read log file; for  species and eggnog level:  "+currentSpe_TaxID +";"+EggNOG_maxLevel
    #     with smtplib.SMTP(smtp_server, port) as server:
    #         server.ehlo()  # Can be omitted
    #         server.starttls(context=context)
    #         server.ehlo()  # Can be omitted
    #         server.login(sender_email, password)
    #         server.sendmail(sender_email, receiver_email, message)




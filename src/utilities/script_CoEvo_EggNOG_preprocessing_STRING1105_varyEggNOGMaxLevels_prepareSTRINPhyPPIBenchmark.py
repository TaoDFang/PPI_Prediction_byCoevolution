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

sys.path.append("../src/utilities")

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
from collect_topCoEvos import get_maxBetValue_dict_pydcaFNAPC_array_npz
from collect_topCoEvos import get_maxBetValue_dict_MI_apc_allResidues_npz
from collect_topCoEvos import getFrame2Dict_bet


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




    try: 

        parser = argparse.ArgumentParser(description='script_CoEvo_EggNOG_preprocessing_STRING11_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark')
        parser.add_argument('-l','--EggNOG_maxLevel', type=str, help='EggNOG_maxLevel')
        parser.add_argument('-i','--currentSpe_TaxID', type=str, help='currentSpe_TaxID')

        args = parser.parse_args()
        EggNOG_maxLevel=args.EggNOG_maxLevel
        currentSpe_TaxID=args.currentSpe_TaxID


        origSTRING_root_folder="/mnt/mnemo6/tao/STRING_Data_11.5/"
        origProSeqPath=origSTRING_root_folder+"STRINGSequencesBySpecies/"
        origSTRINGBacteriaProSeqPath=origSTRING_root_folder+"STRINGBacteriaSequencesBySpecies/"



        newSTRING_rootFolder="/mnt/mnemo6/tao/PPI_Coevolution/STRING_data_11.5/"
        currentSpeProSeqPath=newSTRING_rootFolder+currentSpe_TaxID+"/"
        if not os.path.exists(currentSpeProSeqPath):
            os.makedirs(currentSpeProSeqPath)

        currentSpe_fastaData=currentSpeProSeqPath+currentSpe_TaxID+".fa"
        currentSpeProSeqPath_ByProteins=newSTRING_rootFolder+currentSpe_TaxID+"ByProteins/"


        currentSpeMiddleDataPath=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_MiddleData/"
        if not os.path.exists(currentSpeMiddleDataPath):
            os.makedirs(currentSpeMiddleDataPath)


        # In[ ]:


        origSTRINGBacteriaFiles=glob.glob(origSTRINGBacteriaProSeqPath+"*.fa")
        print("len(origSTRINGBacteriaFiles):",len(origSTRINGBacteriaFiles))


        # In[ ]:





        # # first need to preapre Scerevisiae protein seq data (fasta format )

        # In[ ]:



        #code apdate from http://localhost:8206/lab/workspaces/auto-y/tree/code/STRING_TAO/PPI_Coevolution/Update_STRINGRBH_Scripts/BLAST_CoEvo_preprocessing_qiden.ipynb

        #cp protein seq of current species  to a new folder and separated them by proteins for later use 

        if not os.path.exists(currentSpe_fastaData):
            cmd = ["cp", origProSeqPath+currentSpe_TaxID+".fa", currentSpeProSeqPath] 
            skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            p = subprocess.run(cmd, universal_newlines=True, **skw)
            p.check_returncode()


        # In[ ]:


        # create .fai inndex file for samtool faxid later 
        # to prevent error cause by fun_newSingleMSA_EggNOG_OrthologousGroup_faidx
        # check search fun name in google doc "Things need to do before next meeting 2.0"
        if not os.path.exists(currentSpe_fastaData+".fai"):
            cmd = ["samtools", "faidx",currentSpe_fastaData] # check "python/conda  environment in Jupyterlab" in "PhD环境配置"
            skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            p = subprocess.run(cmd, universal_newlines=True, **skw)
            p.check_returncode()


        # # section ParseEcoliFastaByProtein, need in phmmer section to 

        # In[ ]:


        #%%time
        currentSpe_protein_info_filename=newSTRING_rootFolder+currentSpe_TaxID+".protein.info.v11.5.txt.gz"
        if not os.path.exists(currentSpe_protein_info_filename):
            cmd = [
                "wget", 
                "https://stringdb-static.org/download/protein.info.v11.5/"+currentSpe_TaxID+".protein.info.v11.5.txt.gz",
                '-P', newSTRING_rootFolder

            ]
            skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p = subprocess.run(cmd, universal_newlines=True, **skw)
            p.check_returncode()

        currentSpe_protein_info_frame=pd.read_csv(currentSpe_protein_info_filename,skiprows=1,
                                                  header=None,index_col=None,sep="\t")
        print("currentSpe_protein_info_frame.shape",currentSpe_protein_info_frame.shape)




        if not os.path.exists(currentSpeProSeqPath_ByProteins):
            os.makedirs(currentSpeProSeqPath_ByProteins)
            # this steps some times has problems of get all single proteins , not sur why so far
            ParseCurSpeFastaByProteins_STRNG1105(currentSpe_fastaData,currentSpeProSeqPath_ByProteins)


        print("currentSpeProSeqPath_ByProteins:",currentSpeProSeqPath_ByProteins)
        pro_num=glob.glob(currentSpeProSeqPath_ByProteins+"*")
        print("len(pro_num):",len(pro_num))

        # ParseCurSpeFastaByProteins_STRNG1105 function sometimes dont correctly sperate all single proteins 
        # here to make they do  , otherwise throw error 
        assert currentSpe_protein_info_frame.shape[0]== len(pro_num)

        # CPU times: user 267 ms, sys: 668 ms, total: 935 ms
        # Wall time: 32.3 s


        # In[ ]:


        #currentSpeProSeqPath_ByProteins


        # In[ ]:





        # # remove rudadant proteins for later use 

        # In[ ]:



        #add a new section to remove reduadant  proteins by  aligning  these proteins to each other and 
        #removed redundant ones: only the longer sequence was kept 
        #if two sequences were over 95% identical and the alignment covered 90% of the shorter sequence

        #this can be done by  following :

        #“CreatDBForEcoli” section , run this by copy it to bash shell 


        currentSpeProSeqPath_DB=newSTRING_rootFolder+currentSpe_TaxID+"DB/"+currentSpe_TaxID

        if not os.path.exists(currentSpeProSeqPath_DB):
            os.makedirs(currentSpeProSeqPath_DB)
            cmd = ["/mnt/mnemo5/tao/BeeiveProgram/ncbi-blast-2.10.0+/bin/makeblastdb",
                   "-in", currentSpe_fastaData,
                   "-dbtype","prot",
                   "-out",currentSpeProSeqPath_DB,
                   "-parse_seqids"

                  ] 
            skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            p = subprocess.run(cmd, universal_newlines=True, **skw)
            p.check_returncode()


        currentSpe_withinBlastPath=newSTRING_rootFolder+currentSpe_TaxID+"withinBlast/"

        if not os.path.exists(currentSpe_withinBlastPath):
            os.makedirs(currentSpe_withinBlastPath)
            cmd = ["/mnt/mnemo5/tao/BeeiveProgram/ncbi-blast-2.10.0+/bin/blastp",
                   "-num_threads", "1",
                   "-query",currentSpe_fastaData,
                   "-db",currentSpeProSeqPath_DB,
                   "-out",currentSpe_withinBlastPath+"all2all.txt",
                   "-evalue","1e-6",
                   "-outfmt",'7 qseqid qaccver  qlen sseqid saccver slen qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovs qcovhsp',

                  ] 
            skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            p = subprocess.run(cmd, universal_newlines=True, **skw)
            p.check_returncode()





        # In[ ]:


        if not os.path.exists(newSTRING_rootFolder+currentSpe_TaxID+"_redundant_proteins.csv"):
            withinBlast_result=pd.read_csv(currentSpe_withinBlastPath+"all2all.txt",comment='#',header=None,sep="\t")
            print(withinBlast_result.shape)
            #remove blast reuslt of itself
            withinBlast_result=withinBlast_result.loc[withinBlast_result.iloc[:,0]!=withinBlast_result.iloc[:,3],:]

            print(withinBlast_result.shape)


            withinBlast_result_uniPros=list(set(withinBlast_result.iloc[:,0].values.tolist()))
            print("len(withinBlast_result_uniPros):",len(withinBlast_result_uniPros))
            print(withinBlast_result_uniPros[0:3])


            redudant_tupleList=list() #and the alignment covered 90% of the shorter sequence

            for protein_name  in withinBlast_result_uniPros:
                EcoliwithinBlast_result=withinBlast_result.loc[withinBlast_result.iloc[:,0]==protein_name,:]

                EcoliwithinBlast_result=EcoliwithinBlast_result.loc[EcoliwithinBlast_result.iloc[:,14]>95,:]
                #print(EcoliwithinBlast_result.shape)
                if EcoliwithinBlast_result.shape[0]>0:
                    #print(file_name)
                    for idx in EcoliwithinBlast_result.index:
                        qid,qlen,sid,slen,alignlen=EcoliwithinBlast_result.loc[idx,[0,2,3,5,13]]
                        if qlen<slen:
                            if alignlen>(0.9*qlen):
                                #print(qid)
                                redudant_tupleList.append((sid,qid))

                        else:
                            if alignlen>(0.9*slen):
                                #print(sid)
                                redudant_tupleList.append((qid,sid))


            print("len(redudant_tupleList):",len(redudant_tupleList))
            redudant_tupleList[0:3]




            unique_redudant_tupleList=set(redudant_tupleList)
            print("len(unique_redudant_tupleList):",len(unique_redudant_tupleList))

            redudant_short_pros=[temp[1] for temp in unique_redudant_tupleList] # notice here temp[1] are already shorter redunat protein
            print("len(redudant_short_pros)",len(redudant_short_pros))

            unique_redudant_short_pros=set(redudant_short_pros)
            print("len(unique_redudant_short_pros)",len(unique_redudant_short_pros))
            list(unique_redudant_short_pros)[0:3]



            unique_redudant_short_pros_frame=pd.DataFrame(unique_redudant_short_pros,columns=["pro_name"])
            unique_redudant_short_pros_frame.to_csv(newSTRING_rootFolder+currentSpe_TaxID+"_redundant_proteins.csv",
                                                    header=True,index=None)


        # In[ ]:





        # In[ ]:





        # # then   preprocess eggNOG othologous group , to make sure for each orthologous group  , 
        # only one protein fro one speices 

        # In[ ]:


        #%%time
        # use following code adapated from  "http://localhost:8206/lab/workspaces/auto-I/tree/code/MNF/notebooks/Scerevisiae/choose_orthologs_Scerevisiae.py"
        # now it only works for bacteria 





        currentSpe_currentMaxLevel_orthologs=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_orthologs/"

        if not os.path.exists(currentSpe_currentMaxLevel_orthologs):
            os.makedirs(currentSpe_currentMaxLevel_orthologs)

        #here have to add log file 


            cmd = [
                "python", 
                "/mnt/mnemo5/tao/code/MNF/src/tao_utilities/choose_orthologs_STRING11.05.py",
                '-o',currentSpe_currentMaxLevel_orthologs,
                '-i', currentSpe_fastaData,
                "-m",EggNOG_maxLevel

            ]



            #skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            log_name = os.path.join(currentSpe_currentMaxLevel_orthologs, currentSpe_TaxID + '.log')
            with open(log_name, 'w') as logfp:
                skw = dict(stdout=logfp, stderr=logfp)
                p = subprocess.run(cmd, universal_newlines=True, **skw)
            p.check_returncode()

        # CPU times: user 14.6 ms, sys: 2.01 ms, total: 16.6 ms
        # Wall time: 29min 29s


        # In[ ]:





        # In[ ]:





        # In[ ]:


        # get new othologous group
        generate_ogs_files=glob.glob(currentSpe_currentMaxLevel_orthologs+"*.orthologs")
        print("len(generate_ogs_files):",len(generate_ogs_files))
        generate_ogs_files[0]


        # In[ ]:





        # In[ ]:


        #%%time 
        # check how many of them orthologs files are not empty  or morn then a threshold 
        ortholog_hist=list()
        for test_file in generate_ogs_files:
            test_frame=pd.read_csv(test_file,
                                                 header=None,index_col=None,sep="\t",dtype={0: str,1:str})
            ortholog_hist.append(test_frame.shape[0])


        plt.hist(x=ortholog_hist, bins=40)
        plt.show()


        # In[ ]:





        # In[ ]:


        #%%time 
        # here filter by number of orthologs in single MSA 
        generate_ogs_ArgForCombinedRBH=[(f,len(origSTRINGBacteriaFiles) )for f in generate_ogs_files]
        pool = mp.Pool(30)
        combinedRBH_results = pool.map(combinedRBH_STRING1105, generate_ogs_ArgForCombinedRBH)
        pool.close()





        # In[ ]:


        #%%time 
        print("len(combinedRBH_results):",len(combinedRBH_results))
        combinedRBH_results=[c for c in combinedRBH_results if c is not None]
        print("len(combinedRBH_results):",len(combinedRBH_results))
        combinedRBH_results_flatten=pd.concat(combinedRBH_results)
        print("combinedRBH_results_flatten:",len(combinedRBH_results_flatten))
        combinedRBH_results_flatten.head()


        # In[ ]:





        # In[ ]:





        # In[ ]:





        # #  remove reduant proteins by checking seq similarities 

        # In[ ]:




        #code apdate from http://localhost:8206/lab/workspaces/auto-y/tree/code/STRING_TAO/PPI_Coevolution/Update_STRINGRBH_Scripts/BLAST_CoEvo_preprocessing_qiden.ipynb


        # In[ ]:


        currentSpe_redundant_pros=pd.read_csv(newSTRING_rootFolder+currentSpe_TaxID+"_redundant_proteins.csv",
                                         header=0,index_col=None)
        print("currentSpe_redundant_pros.shape:",currentSpe_redundant_pros.shape)
        currentSpe_redundant_pros_list=list(currentSpe_redundant_pros['pro_name'])
        currentSpe_redundant_pros_list[0:3]


        # In[ ]:




        NonRedundantIdx=[False if pro in currentSpe_redundant_pros_list else True for pro in combinedRBH_results_flatten.iloc[:,0].tolist()]


        combinedRBH_results_frame=combinedRBH_results_flatten.rename(columns={0:"query_protein",1:"RBH_hit",})
        print("combinedRBH_results_frame.shape:",combinedRBH_results_frame.shape)
        combinedRBH_results_frame_NonRedundantIdx=combinedRBH_results_frame.loc[NonRedundantIdx,:]
        print("combinedRBH_results_frame_NonRedundantIdx.shape:",combinedRBH_results_frame_NonRedundantIdx.shape)
        combinedRBH_results_frame_NonRedundantIdx.head()


        # In[ ]:


        Unique_currentSpe_Pros=list(set(combinedRBH_results_frame_NonRedundantIdx['query_protein']))
        #Unique_Ecoli_Pros=list(set(Ecoli_within.iloc[:,0]))

        print("len(Unique_currentSpe_Pros),Unique_currentSpe_Pros[0:3]",len(Unique_currentSpe_Pros),Unique_currentSpe_Pros[0:3]) #check this value


        # In[ ]:





        # In[ ]:





        # In[ ]:


        #%%time 

        def getCurrentSpeRBHdict(pro):
            #print(pro)
            RBH_HitsResults_Otholog_dic=defaultdict(list)

            currentSpe_RBH_frame=combinedRBH_results_frame_NonRedundantIdx.loc[combinedRBH_results_frame_NonRedundantIdx['query_protein']==pro,:]

            currentSpe_string_proName=list(set(currentSpe_RBH_frame['RBH_hit']))
            currentSpe_string_speId=[temp.split(".")[0] for temp in currentSpe_string_proName]



            for  idx, forward_hit in enumerate(currentSpe_string_proName):
                if currentSpe_TaxID not in forward_hit: # here add this as damains's file alread contain ecoli protein  add the begining 
                    forward_hit_SpeId=currentSpe_string_speId[idx]
                    RBH_HitsResults_Otholog_dic[pro].append(forward_hit_SpeId+"|"+forward_hit)

            if len(RBH_HitsResults_Otholog_dic.keys())==0:
                print(pro)
            else:
                return(RBH_HitsResults_Otholog_dic)




        newsingleMSA_RBH_OrthologousGroup_fileName=currentSpeMiddleDataPath+"newsingleMSA_RBH_OrthologousGroup.csv"

        if not os.path.exists(newsingleMSA_RBH_OrthologousGroup_fileName):
            pool = mp.Pool(30) # 
            results = pool.map(getCurrentSpeRBHdict, Unique_currentSpe_Pros) # here many damian generalted files are emply 
            pool.close() 

            results=[temp for temp in results if temp is not None]
            print("len(results):",len(results))


            # only run once as its append mode 
            for RBH_HitsResults_Otholog_dic_temp in results:
                with open(newsingleMSA_RBH_OrthologousGroup_fileName,"a") as file:
                    for key in RBH_HitsResults_Otholog_dic_temp.keys():
                        file.write("%s,%s\n" % (key, ",".join(RBH_HitsResults_Otholog_dic_temp[key]))) 


        # CPU times: user 251 ms, sys: 321 ms, total: 573 ms
        # Wall time: 20.3 s


        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # # for each OG group, find their sequence and save them in to one fasta file ,

        # In[ ]:


        #%%time
        #for each OG group, find their sequence and save them in to one fasta file ,
        #first seuquence has to be query/ecoli/Scerevisiae/  protein for downstream filtering 
        # this step is quite slow so before its beeive were used but now we want whole script can run in same place 





        currentSpe_OrthologousGroup_Fa_path=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_EggNOG_OrthologousGroup_Fa/"

        currentSpe_OrthologousGroup_Fa_logpath=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_EggNOG_OrthologousGroup_Fa_log/"

        if not os.path.exists(currentSpe_OrthologousGroup_Fa_path):
            os.makedirs(currentSpe_OrthologousGroup_Fa_path) 
            os.makedirs(currentSpe_OrthologousGroup_Fa_logpath) 



            newsingleMSA_RBH_OrthologousGroup_frame=pd.read_csv(newsingleMSA_RBH_OrthologousGroup_fileName,
                                                                header=None,index_col=None, sep=" ")
            print(newsingleMSA_RBH_OrthologousGroup_frame.shape)
            newsingleMSA_RBH_OrthologousGroup_OGidx=[idx+1 for idx in range(newsingleMSA_RBH_OrthologousGroup_frame.shape[0])]

            newsingleMSA_RBH_OrthologousGroup_OGidx_ArgForFaidx=[(currentSpe_TaxID,OG_idx,newSTRING_rootFolder,origSTRINGBacteriaProSeqPath,currentSpe_OrthologousGroup_Fa_path,currentSpe_OrthologousGroup_Fa_logpath,newsingleMSA_RBH_OrthologousGroup_fileName) for OG_idx in newsingleMSA_RBH_OrthologousGroup_OGidx]
            pool=mp.Pool(30)  # here when we set it as 50, there are more i.o wait time from top command ; 
            pool.map(fun_newSingleMSA_EggNOG_OrthologousGroup_faidx,newsingleMSA_RBH_OrthologousGroup_OGidx_ArgForFaidx)
            pool.close() 


        # CPU times: user 6.74 s, sys: 43.2 s, total: 49.9 s
        # Wall time: 33min 25s


        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # # “SeedAlignment” section 
        # See suplemenaty files of Qian's paper : 
        # Phmmer choose most similar sequce
        # and hmmalig to do multiple alignment ?
        # http://www.csb.yale.edu/userguides/seq/hmmer/docs/node7.html 
        # and read Pfam paper

        # In[ ]:


        #%%time 

            
            
        currentSpe_OrthologousGroup_Fa_files=glob.glob(currentSpe_OrthologousGroup_Fa_path+"*")
        print("len(currentSpe_OrthologousGroup_Fa_files):",len(currentSpe_OrthologousGroup_Fa_files))
        currentSpe_OrthologousGroup_Fa_idxs=[i+1 for i in range(len(currentSpe_OrthologousGroup_Fa_files))]

        currentSpe_phmmer_outPath=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_phmmer_results_pident/"
        currentSpe_phmmer_logPath=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_phmmer_results_pident_log/"

        if not os.path.exists(currentSpe_phmmer_outPath):
            os.makedirs(currentSpe_phmmer_outPath) 
            os.makedirs(currentSpe_phmmer_logPath) 
            currentSpe_OrthologousGroup_Fa_idxs_ArgFroPhmmer=[(currentSpe_TaxID,idx,currentSpe_OrthologousGroup_Fa_path,currentSpeProSeqPath_ByProteins,currentSpe_phmmer_outPath,currentSpe_phmmer_logPath) for idx in currentSpe_OrthologousGroup_Fa_idxs]
            pool=mp.Pool(50) 
            pool.map(fun_newSingleMSA_EggNOG_OrthologousGroup_phmmer,currentSpe_OrthologousGroup_Fa_idxs_ArgFroPhmmer)
            pool.close() 


        # CPU times: user 1.51 s, sys: 35 s, total: 36.5 s
        # Wall time: 6min 54s



        # In[ ]:





        # # then filter phmmer resutls using next two  blocks
        # 
        # 

        # In[ ]:


        #%%time

        def getSpe2pro2seq_dict(record):
            origSTRINGBacteriaProSeqPath,SpeName=record
            Spe_index=SeqIO.index(origSTRINGBacteriaProSeqPath+SpeName+".fa","fasta")
            Spe_index_dict=dict(Spe_index)
            Spe_index.close()
            return((SpeName,Spe_index_dict))


        # keep this funcion heree as it need very larget global variable 
        # and using subprocess seem not very fast 
        def newSingleMSA_filterPhmmer_EggNOG_OrthologousGroup_faidx(file_name):
            try:
                reTemp=re.compile("([0-9]*)\.*")
                #currentSpe_phmmer_outPath,currentSpe_phmmer_OrthologousGroup_path,file_name,STRING_id2newname_dict,species_indexDict_dict=record
                init_phmmer_result=pd.read_csv(currentSpe_phmmer_outPath+file_name, 
                                    comment="#", delim_whitespace=True,header=None,index_col=None,engine="python")
                init_phmmer_result['qcov']=(init_phmmer_result.iloc[:,16]-init_phmmer_result.iloc[:,15]+1)/init_phmmer_result.iloc[:,5] # query coverage
                init_OGLen=init_phmmer_result.shape[0]
                # targetname accession   tlen queryname       0-3
                #accession   qlen   E-value  score  bias   4-8
                # #  of  c-Evalue  i-Evalue  score  bias  9-14
                # from    to  from    to    #15-18
                #from    to  acc  description oftarget #19-22
                # checck Tabular output formats in hmmer documentation 
                # here acc is not identitiy .but a similar one
                phmmer_result=init_phmmer_result.loc[init_phmmer_result.iloc[:,21]>0.55,:].copy() # acc > 0.55 acc >0.4
                phmmer_result=phmmer_result.loc[phmmer_result['qcov']>0.8,:] # query coverage > 0.8 >0.65
                if phmmer_result.shape[0]>2500 or phmmer_result.shape[0]>init_OGLen*0.25:
                    pass
                else:
                    #print("second layer")
                    phmmer_result=init_phmmer_result.loc[init_phmmer_result.iloc[:,21]>0.4,:].copy() # acc > 0.55 acc >0.4
                    phmmer_result=phmmer_result.loc[phmmer_result['qcov']>0.65,:] # query coverage > 0.8 >0.65
                    if phmmer_result.shape[0]>2500 or phmmer_result.shape[0]>init_OGLen*0.25:
                        pass
                    else:
                        #print("third layer")
                        phmmer_result=init_phmmer_result.loc[init_phmmer_result.iloc[:,21]>0.25,:].copy() # acc > 0.55 acc >0.4
                        phmmer_result=phmmer_result.loc[phmmer_result['qcov']>0.5,:] # query coverage > 0.8 >0.65

                phmmer_result=phmmer_result.loc[phmmer_result.iloc[:,0]!=file_name[0:-10],:] 
            except:
                print("empty file : "+file_name)
                #print(file_name)



            #if (phmmer_result.shape[0]>0):
            if (phmmer_result.shape[0]>2500 or phmmer_result.shape[0]>init_OGLen*0.25):
                currentSpe_protein=phmmer_result.iloc[0,3]
                ecoli_dict=SeqIO.index(currentSpeProSeqPath_ByProteins+currentSpe_protein+".fa","fasta")
                ecoli_str=str(ecoli_dict[currentSpe_protein].seq)
                output_str=">"+currentSpe_protein+"\n"+ecoli_str+"\n"

                phmmer_ecoli=phmmer_result.iloc[0,3]
                #print(phmmer_ecoli)
                #???!!! here there are duplicete target proteins as it has  more then two domains aligned 
                # but here we only use unique targert prorteins name and check if duplicted later  
                phmmer_OG_Pros=list(set(phmmer_result.iloc[:,0]))


            # Here is a little tricky
            #which is suppose to slover the problem mention in page 21 of suplemntary files of Qians paperif 
            #First we need to consider alignment(target) overlap in query,  if so, only take on with high identity
            #Second we need to condiser alignent(target) overlap in target, if so. useing intersection of them.
            #or orignal sequence of them. if so. the target sequececn len will increase!!
            #how about only conside the situatin when non overlap in both query and target sequence ?
                for hit_seq in phmmer_OG_Pros: # 
                    temp_phmmer_result=phmmer_result.loc[phmmer_result.iloc[:,0]==hit_seq,:].copy()
                    if temp_phmmer_result.shape[0]>1:
                        temp_phmmer_result.sort_values(by=[15],ascending=True,inplace=True)# make sure segment in query in correct order
                        current_index=0
                        keepedIdx=list()
                        for i in range(1,temp_phmmer_result.shape[0]): # get non overlaped  and order segment in query 
                    #         if i>range(temp_phmmer_result.shape[0]):
                    #             break
                    #        else:
                            if temp_phmmer_result.iloc[current_index,15]<=temp_phmmer_result.iloc[i,15]<=temp_phmmer_result.iloc[current_index,16]:# overlap 
                                #print("Im here")
                                if temp_phmmer_result.iloc[current_index,21]<temp_phmmer_result.iloc[i,21]:
                                    #print("Im there")
                                    current_index=i
                                if i==(temp_phmmer_result.shape[0]-1) or (temp_phmmer_result.iloc[i+1,15]>temp_phmmer_result.iloc[current_index,16]):
                                    keepedIdx.append(current_index)
                            elif temp_phmmer_result.iloc[i,15]>temp_phmmer_result.iloc[current_index,16]:   # none overlap and next segment is after current segemtn 
                          #      keepedIdx.append(current_index) # thi can be omited 
                                keepedIdx.append(i)
                                current_index=i

                        keepedIdx=sorted(list(set(keepedIdx)))
                        temp_temp_phmmer_result=temp_phmmer_result.iloc[keepedIdx,:].copy()


                        current_index2=0
                        keepedIdx2=[0]
                        for i in range(1,temp_temp_phmmer_result.shape[0]): # get non overlaped and order segments in target
            #                 if i>range(temp_temp_phmmer_result.shape[0]):
            #                     break
            #                 else:
                            if temp_temp_phmmer_result.iloc[current_index,18]<temp_temp_phmmer_result.iloc[i,17]:
                                current_index2=i
                                keepedIdx2.append(current_index2)
                        final_temp_phmmer_result=temp_temp_phmmer_result.iloc[keepedIdx2,:]
                    else: 
                        final_temp_phmmer_result=temp_phmmer_result

                    hit_str=str()
                    if final_temp_phmmer_result.shape[0]>1:
                        print("im non overlaped in both query and target: "+file_name)
                        #print(file_name)
                    for j in range(final_temp_phmmer_result.shape[0]):
                        hit_seq,hit_start,hit_end=final_temp_phmmer_result.iloc[j,[0,17,18]]
                        hit_id=reTemp.match(hit_seq).group(1)
                        #hit_SpeName=STRING_id2newname_dict[int(hit_id)]

                        #hit_SpeName=STRING_id2newname.loc[int(hit_id),'name']
                        #print(hit_SpeName)
                        # originall Ecoli is not added as first line so this folder is oaky,
                        # but with now pipeline , i need to move ecoli species back to /mnt/mnemo5/tao/PPI_Coevolution/Uniprot_Ecoli/Updated_STRINGRBH_pipeline/STRINGBacteriaSequencesBySpecies
                        #(base) tao@gaia:~/STRING$ cp STRINGSequencesBySpecies/Escherichia_coli_str._K-12_substr._MG1655.fa  STRINGBacteriaSequencesBySpecies/
                        record_dict = species_indexDict_dict[hit_id]
                        hit_str=hit_str+str(record_dict[hit_seq].seq)[(hit_start-1):(hit_end-1)]

                    output_str=output_str+">"+hit_seq+"\n"+hit_str+"\n"
                        #phmmer_Otholog_dic[phmmer_ecoli].append(hit_SpeName+"|"+hit_seq)
            #         with open("/mnt/mnemo5/tao/PPI_Coevolution/Uniprot_Ecoli/Updated_STRINGRBH_pipeline/phmmer_OrthologousGroup_Fa_pident/"+currentSpe_protein+".fa","w") as out_str_file:
            #             out_str_file.write(output_str) # whne use with open, print_time not working. because not imediately close??
                out_str_file=open(currentSpe_phmmer_OrthologousGroup_path+currentSpe_protein+".fa","w")
                out_str_file.write(output_str)
                out_str_file.close()

            else:
                print("not enough seq in phmmer seq alignment ",file_name+":",phmmer_result.shape)


        currentSpe_phmmer_OrthologousGroup_path=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_phmmer_OrthologousGroup_Fa_pident/"

        if not os.path.exists(currentSpe_phmmer_OrthologousGroup_path):
            os.makedirs(currentSpe_phmmer_OrthologousGroup_path)


            species_fa=glob.glob(origSTRINGBacteriaProSeqPath+"*.fa")
            species_names=[os.path.basename(f) for f in species_fa]
            species_names=[f[:-3] for f in species_names]
            print("len(species_names):",len(species_names))




            species_names_ArgForGetSpe2pro2seq=[(origSTRINGBacteriaProSeqPath,SpeName) for SpeName in species_names]
            pool = mp.Pool(10)
            getSpe2pro2seq_dict_results=pool.map(getSpe2pro2seq_dict, species_names_ArgForGetSpe2pro2seq)
            pool.close()

            species_indexDict_dict=dict(getSpe2pro2seq_dict_results)


            # CPU times: user 4min 58s, sys: 30.2 s, total: 5min 28s
            # Wall time: 5min 36s





            #run once 
            # now for one species around 3-4 minutes
            # here many operation involde with find one particular protei seq from on species, 
            # better to read them in one run 
            # after read all spefice fasta file, one run takes  3.81 seconds 


            currentSpe_phmmer_files = [f for f in listdir(currentSpe_phmmer_outPath) if isfile(join(currentSpe_phmmer_outPath, f))]
            currentSpe_phmmer_files=[f for f in currentSpe_phmmer_files if "domtblout" in f]

            print("len(currentSpe_phmmer_files):",len(currentSpe_phmmer_files))

            #currentSpe_phmmer_files_ArgForfilterPhmmer=[(currentSpe_phmmer_outPath,currentSpe_phmmer_OrthologousGroup_path,file_name,STRING_id2newname_dict,species_indexDict_dict) for file_name in currentSpe_phmmer_files]
            pool = mp.Pool(10) # here dont use many process here as following functin require large memory 
            pool.map(newSingleMSA_filterPhmmer_EggNOG_OrthologousGroup_faidx, currentSpe_phmmer_files)
            pool.close()

        # CPU times: user 5min 45s, sys: 37.8 s, total: 6min 23s
        # Wall time: 13min 46s


        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # # “ProfileHMM” section

        # In[ ]:


        # before to make HMM profile from seed alignmt. we first need to perfomacne multiple sequence alignment 


        # In[ ]:


        #%%time 

        currentSpe_ClustoMSA_path=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_phmmer_OrthologousGroup_Fa_pident_ClustoMSA/"

        if not os.path.exists(currentSpe_ClustoMSA_path):
            os.makedirs(currentSpe_ClustoMSA_path)
            phmmer_OrthologousGroup_names=glob.glob(currentSpe_phmmer_OrthologousGroup_path+"*.fa")
            phmmer_OrthologousGroup_names=[os.path.basename(f) for f in phmmer_OrthologousGroup_names]
            phmmer_OrthologousGroup_names=[f[:-3] for f in phmmer_OrthologousGroup_names]

            print("len(phmmer_OrthologousGroup_names),phmmer_OrthologousGroup_names[0]:",len(phmmer_OrthologousGroup_names),phmmer_OrthologousGroup_names[0])


            phmmer_OrthologousGroup_names_ArgFroClustoMSA=[(name,currentSpe_phmmer_OrthologousGroup_path,currentSpe_ClustoMSA_path) for name in phmmer_OrthologousGroup_names]
            pool = mp.Pool(50)
            pool.map(ClustoMSA, phmmer_OrthologousGroup_names_ArgFroClustoMSA)
            pool.close()



        # CPU times: user 34.3 s, sys: 1min 34s, total: 2min 9s
        # Wall time: 1h 11min 7s


        # In[ ]:





        # In[ ]:


        # check columsn/aligned postions in ClustoO MSA file 
        Clusto=glob.glob(currentSpe_ClustoMSA_path+"*")
        print("len(Clusto):",len(Clusto))
        print(Clusto[0])


        # In[ ]:


        #%%time 

        #Better draw a histogram!!! 
        def getClustoMSAlen(file_name):
            Clusto_msa=AlignIO.read(file_name,"stockholm")
            return(Clusto_msa.get_alignment_length())


        pool = mp.Pool(50)
        Clusto_msa_hist=pool.map(getClustoMSAlen, Clusto)
        pool.close()


        plt.hist(Clusto_msa_hist,bins=50)
        plt.show()

        # CPU times: user 4.55 ms, sys: 44.4 s, total: 44.4 s
        # Wall time: 33.8 s


        # In[ ]:


        min(Clusto_msa_hist),max(Clusto_msa_hist)


        # In[ ]:





        # # now build hmm profile use seed alignment  using next block
        # 
        # 

        # In[ ]:


        #%%time 

        currentSpe_hmm_profiles_path=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_hmm_profiles_pident/"

        if not os.path.exists(currentSpe_hmm_profiles_path):
            os.makedirs(currentSpe_hmm_profiles_path)

            ClustoMSA_files=glob.glob(currentSpe_ClustoMSA_path+"*")
            ClustoMSA_files=[os.path.basename(f) for f in ClustoMSA_files]
            print("len(ClustoMSA_files):",len(ClustoMSA_files),ClustoMSA_files[0])

            ClustoMSA_files_ArgForhmmbuild=[(name,currentSpe_hmm_profiles_path,currentSpe_hmm_profiles_path,currentSpe_ClustoMSA_path) for name in ClustoMSA_files]
            pool = mp.Pool(50)
            pool.map(hmmbuild, ClustoMSA_files_ArgForhmmbuild)
            pool.close()

        # CPU times: user 370 ms, sys: 44.7 s, total: 45.1 s
        # Wall time: 2min 38s



        # In[ ]:





        # In[ ]:





        # In[ ]:





        # # “HMMAllAlignment” using next block 
        # 
        # 

        # In[ ]:


        #%%time 




        currentSpe_hmm_align_path=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_hmmalign_pident_results/"

        if not os.path.exists(currentSpe_hmm_align_path):
            os.makedirs(currentSpe_hmm_align_path)

            hmm_files=glob.glob(currentSpe_hmm_profiles_path+"*.hmm")
            hmm_files=[os.path.basename(f) for f in hmm_files]
            hmm_files=[f[:-4] for f in hmm_files]

            print("len(hmm_files):",len(hmm_files),hmm_files[0])


            hmm_files_ArgForhmmalign=[(name,currentSpe_hmm_align_path,currentSpe_hmm_profiles_path,currentSpe_OrthologousGroup_Fa_path) for name in hmm_files]
            pool = mp.Pool(50)
            pool.map(hmmalign, hmm_files_ArgForhmmalign)
            pool.close()

        # CPU times: user 1.2 s, sys: 1.37 s, total: 2.57 s
        # Wall time: 9min 3s


        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # # Here we need to further remove gaps from single MSA and keep track their position in original seq 
        # here we first remove gaps in ref seq / first protein 
        # then remove seq with many gaps 
        # then remove columns with many gaps 
        # 

        # In[ ]:


        #%%time 


        currentSpe_msa_removeGaps_path=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_hmmalign_removeGaps_keepGapPos/"
        currentSpe_msa_trackGapsPos_path=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_hmmalign_removeGaps_trackGapPos/"



        if not os.path.exists(currentSpe_msa_trackGapsPos_path):
            os.makedirs(currentSpe_msa_trackGapsPos_path)

        if not os.path.exists(currentSpe_msa_removeGaps_path):
            os.makedirs(currentSpe_msa_removeGaps_path)

        #  seem sthis can killl kernel 
        #     import warnings # https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing
        #     warnings.filterwarnings('error')


            hmmalign_sto_files=glob.glob(currentSpe_hmm_align_path+"/*.sto")
            print("len(hmmalign_sto_files):",len(hmmalign_sto_files))

            #all_idx=list(range(len(hmmalign_sto_files)))
            hmmalign_sto_files_ArgForremoeGaps=[(file_name,currentSpe_ClustoMSA_path,currentSpe_msa_trackGapsPos_path,currentSpe_msa_removeGaps_path) for file_name in hmmalign_sto_files]
            pool = mp.Pool(40)
            pool.map(removeGapsANDtrackAAPosOfSeqInMSA, hmmalign_sto_files_ArgForremoeGaps)
            pool.close()


        # CPU times: user 3.73 s, sys: 3.18 s, total: 6.9 s
        # Wall time: 17min 22s


        # In[ ]:





        # In[ ]:





        # # get final single MSA len
        # 

        # In[ ]:





        # In[ ]:


        final_MSAs=glob.glob(currentSpe_msa_removeGaps_path+"*")
        print("len(final_MSAs):",len(final_MSAs))
        print(final_MSAs[0])

        final_MSA_pros=[os.path.basename(f) for f in final_MSAs]
        final_MSA_pros=[f[:-3] for f in final_MSA_pros]
        print(final_MSA_pros[0])


        # In[ ]:


        #%%time
        # get len and Nf90 of all single MSA 
        if not os.path.exists(currentSpeMiddleDataPath+'fasta_protein_lens_dict.pickle'):

            fasta_protein_lens=dict()
            fasta_protein_Nf90s=dict()

            for pro in final_MSA_pros:
                alig = AlignIO.read(currentSpe_msa_removeGaps_path+pro+".fa", "fasta")
                fasta_protein_lens[pro]=len(alig[0])
                fasta_protein_Nf90s[pro]=len(alig)/math.sqrt(len(alig[0]))

            with open(currentSpeMiddleDataPath+'fasta_protein_lens_dict.pickle', 'wb') as handle:
                pickle.dump(fasta_protein_lens, handle)

            with open(currentSpeMiddleDataPath+'fasta_protein_Nf90s_dict.pickle', 'wb') as handle:
                pickle.dump(fasta_protein_Nf90s, handle)


        # In[ ]:





        # In[ ]:





        # # -----------------prepare benchmark dataset -----------------------------
        # 
        # code apate from http://localhost:8306/lab/tree/code/MNF/notebooks/Extend2OtherBacterialSpe/CoEvo_EggNOG_prepareSTRINPhyPPIBenchmark.ipynb
        # 

        # In[ ]:





        # # define some path and read some necesseary data 

        # In[ ]:





        # In[ ]:



        CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING11.5/"


        input_root_folder=CoEvo_data_folder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_eggNOGfilteredData/"

        Benchmark_folder=input_root_folder+"STRINPhyPPI_Benchmark/"


        pairedMSA_unfiltered_folder=input_root_folder+"pair_MSA_unfiltered_PasteAlign/"
        pairedMSA_hhfilter_folder=input_root_folder+"pair_MSA_hhfilter_PasteAlign/"
        pairedMSA_Nf90_folder = input_root_folder+"pair_MSA_Nf90_PasteAlign/"
        pairedMSA_Nf90_csv = input_root_folder+"pair_MSA_Nf90.csv"
        pairedMSA_sameProteinRatio_csv = input_root_folder+"sameProteinRatio.csv"

        DCA_coevolutoin_path=input_root_folder+"coevolutoin_result_DCA/"
        MI_coevolutoin_path=input_root_folder+"coevolutoin_result_MI/"

        if not os.path.exists(input_root_folder):
            os.makedirs(input_root_folder)

        if not os.path.exists(Benchmark_folder):
            os.makedirs(Benchmark_folder)

        if not os.path.exists(DCA_coevolutoin_path):
            os.makedirs(DCA_coevolutoin_path)

        if not os.path.exists(MI_coevolutoin_path):
            os.makedirs(MI_coevolutoin_path)

        if not os.path.exists(pairedMSA_unfiltered_folder):
            os.makedirs(pairedMSA_unfiltered_folder)

        if not os.path.exists(pairedMSA_hhfilter_folder):
            os.makedirs(pairedMSA_hhfilter_folder)

        if not os.path.exists(pairedMSA_Nf90_folder):
            os.makedirs(pairedMSA_Nf90_folder)

        if not os.path.exists(pairedMSA_Nf90_csv):
            with open(pairedMSA_Nf90_csv, 'w') as fp:
                pass

        if not os.path.exists(pairedMSA_sameProteinRatio_csv):
            with open(pairedMSA_sameProteinRatio_csv, 'w') as fp:
                pass




        currentSpe_hmmalign_path=newSTRING_rootFolder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_newSingleMSA_hmmalign_removeGaps_keepGapPos/"




        # In[ ]:


        Nf90_thres=16

        with open(currentSpeMiddleDataPath+'fasta_protein_lens_dict.pickle', 'rb') as handle:
            fasta_protein_lens=pickle.load(handle)

        with open(currentSpeMiddleDataPath+'fasta_protein_Nf90s_dict.pickle', 'rb') as handle:
            fasta_protein_Nf90s=pickle.load(handle)




        currentSpe_final_MSAs=glob.glob(currentSpe_hmmalign_path+"*")
        print("len(currentSpe_final_MSAs):",len(currentSpe_final_MSAs))
        print(currentSpe_final_MSAs[0])

        currentSpe_final_pros=[os.path.basename(f) for f in currentSpe_final_MSAs]
        currentSpe_final_pros=[p[:-3] for p in currentSpe_final_pros]
        print(currentSpe_final_pros[0])

        currentSpe_final_pro_dict=dict([(p,1)for p in currentSpe_final_pros])


        # In[ ]:


        #%%time
        # clean keeg benchmakr datset with string score 

        # this is a reversed file alread :(p1,p2),(p2,p1)
        currentSpe_string_score_filename=newSTRING_rootFolder+currentSpe_TaxID+".protein.links.detailed.v11.5.txt.gz"
        if not os.path.exists(currentSpe_string_score_filename):
            cmd = [
                "wget", 
                "https://stringdb-static.org/download/protein.links.detailed.v11.5/"+currentSpe_TaxID+".protein.links.detailed.v11.5.txt.gz",
                '-P', newSTRING_rootFolder

            ]



            skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p = subprocess.run(cmd, universal_newlines=True, **skw)
            p.check_returncode()


        currentSpe_string_score=pd.read_csv(currentSpe_string_score_filename,
                                       header=0,index_col=None,sep=" ")
        currentSpe_string_score_list=currentSpe_string_score.loc[:,["protein1","protein2","combined_score"]].values.tolist()
        currentSpe_string_score_dict=dict([((p1,p2),s)for p1, p2, s in currentSpe_string_score_list])

        print("(currentSpe_string_score.shape:",currentSpe_string_score.shape)


        # CPU times: user 47.2 s, sys: 24.2 s, total: 1min 11s
        # Wall time: 1min 12s


        # In[ ]:





        # In[ ]:





        # # prepare benchmark dataset 
        # 

        # 1: prepare positive and negative ppi 
        # 2: need to prepair paired MSA ,write one function to process one protein pair 
        # 3: need to calculate DCA or MI data 
        # 4: further combined info mation 

        # ## get postivte and negative ppi  first step 

        # In[ ]:


        STRINGcurrentSpePhyPPI_benchmark_file="/mnt/mnemo6/damian/STRING_derived_v11.5/download_files/protein.physical.links.v11.5/"+currentSpe_TaxID+".protein.physical.links.v11.5.txt.gz"
        STRINGcurrentSpePhyPPI_benchmark=pd.read_csv(STRINGcurrentSpePhyPPI_benchmark_file,sep=" ",
                                       header=0,index_col=None)

        print("STRINGcurrentSpePhyPPI_benchmark.shape:",STRINGcurrentSpePhyPPI_benchmark.shape)

        STRINGcurrentSpePhyPPI_benchmark=STRINGcurrentSpePhyPPI_benchmark.loc[STRINGcurrentSpePhyPPI_benchmark['combined_score']>500,:]
        print("STRINGcurrentSpePhyPPI_benchmark.shape:",STRINGcurrentSpePhyPPI_benchmark.shape)

        STRINGcurrentSpePhyPPI_benchmark=STRINGcurrentSpePhyPPI_benchmark.sort_values(by="combined_score",ascending=False)

        STRINGcurrentSpePhyPPI_benchmark.head(n=3)



        #set(STRINGcurrentSpePhyPPI_benchmark.iloc[:,1])

        # notice here  already contains reversed version ppi 


        currentSpe_STRINGcurrentSpePhyPPI_posPPI=STRINGcurrentSpePhyPPI_benchmark.loc[:,["protein1","protein2"]].values.tolist()



        # here rememver one thing different 
        # if we start to calculate dca score only for ppi in benchmark dataset, we shouldnt use use reversed version of pp
        # its a litter different as we do in Ecoli_KEGG , where we first calcuate dca for all protein pari in this species 
        currentSpe_STRINGcurrentSpePhyPPI_posPPI=[tuple(sorted(pp)) for pp in currentSpe_STRINGcurrentSpePhyPPI_posPPI]
        currentSpe_STRINGcurrentSpePhyPPI_posPPI=list(set(currentSpe_STRINGcurrentSpePhyPPI_posPPI))


        currentSpe_STRINGcurrentSpePhyPPI_posPPI=[tuple(ppi) for ppi in currentSpe_STRINGcurrentSpePhyPPI_posPPI]
        # here need to make sure both protein have alignment data 
        currentSpe_STRINGcurrentSpePhyPPI_posPPI=[(p1,p2)for p1, p2 in currentSpe_STRINGcurrentSpePhyPPI_posPPI if (p1 in currentSpe_final_pro_dict) and (p2 in currentSpe_final_pro_dict)]

        print("len(currentSpe_STRINGcurrentSpePhyPPI_posPPI):",len(currentSpe_STRINGcurrentSpePhyPPI_posPPI))
        print(currentSpe_STRINGcurrentSpePhyPPI_posPPI[0:3])


        currentSpe_STRINGcurrentSpePhyPPI_posPPI=sorted(list(set(currentSpe_STRINGcurrentSpePhyPPI_posPPI)))
        print("len(currentSpe_STRINGcurrentSpePhyPPI_posPPI)",len(currentSpe_STRINGcurrentSpePhyPPI_posPPI))


        currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict=dict([(pp,1) for pp in currentSpe_STRINGcurrentSpePhyPPI_posPPI])
        print("len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict)",len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict))
        print(list(currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict.items())[0:3])




        currentSpe_STRINGcurrentSpePhyPPI_allProteins=[p for key in currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict.keys() for p in key ]

        currentSpe_STRINGcurrentSpePhyPPI_allProteins=list(set(currentSpe_STRINGcurrentSpePhyPPI_allProteins))

        # here is import to make sure  every run of this pipeline producing same benchmark 
        currentSpe_STRINGcurrentSpePhyPPI_allProteins= sorted(set(currentSpe_STRINGcurrentSpePhyPPI_allProteins))

        print(currentSpe_STRINGcurrentSpePhyPPI_allProteins[0:3])

        print("len(currentSpe_STRINGcurrentSpePhyPPI_allProteins):",len(currentSpe_STRINGcurrentSpePhyPPI_allProteins))
        currentSpe_STRINGcurrentSpePhyPPI_allProteins_dict=dict([(p,1) for p in currentSpe_STRINGcurrentSpePhyPPI_allProteins])




        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:


        #%%time
        # constructure final positive  PPI 

        currentSpe_STRINGcurrentSpePhyPPI_posPPI_info=[[p1,p2,fasta_protein_lens[p1],fasta_protein_lens[p2]] for p1,p2 in currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict.keys()]
        print("len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_info)",len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_info))



        # In[ ]:


        # constracut negaive PPI using proteins in positive PPI

        currentSpe_STRINGcurrentSpePhyPPI_negPPI=list(itertools.combinations(currentSpe_STRINGcurrentSpePhyPPI_allProteins, 2))
        print("len(currentSpe_STRINGcurrentSpePhyPPI_negPPI)",len(currentSpe_STRINGcurrentSpePhyPPI_negPPI))
        currentSpe_STRINGcurrentSpePhyPPI_negPPI=[(p1,p2) for p1,p2 in currentSpe_STRINGcurrentSpePhyPPI_negPPI if (((p1,p2) not in currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict) and ((p2,p1) not in currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict)) ]


        print("len(currentSpe_STRINGcurrentSpePhyPPI_negPPI):",len(currentSpe_STRINGcurrentSpePhyPPI_negPPI))

        # here shuffle negative ppi randomly 
        random.Random(10).shuffle(currentSpe_STRINGcurrentSpePhyPPI_negPPI)
        print(currentSpe_STRINGcurrentSpePhyPPI_negPPI[0:3])

        # here random sample from negative ppi 

        # here random.sample cant be fixed , so might cause different neg in different runs. better to limit protein frequency later
        # max_negPPI_num=min(len(currentSpe_STRINGcurrentSpePhyPPI_negPPI),len(currentSpe_STRINGcurrentSpePhyPPI_posPPI)*10)
        # currentSpe_STRINGcurrentSpePhyPPI_negPPI=random.sample(currentSpe_STRINGcurrentSpePhyPPI_negPPI,max_negPPI_num)
        # print(currentSpe_STRINGcurrentSpePhyPPI_negPPI[0:3])

        currentSpe_STRINGcurrentSpePhyPPI_negPPI_info=[[p1,p2,fasta_protein_lens[p1],fasta_protein_lens[p2]] for p1,p2 in currentSpe_STRINGcurrentSpePhyPPI_negPPI]
        print("len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_info):",len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_info))


        # In[ ]:





        # In[ ]:


        #%%time
        # constructure final positive and negigeve PPI using string score 

        #currentSpe_STRINGcurrentSpePhyPPI_negPPI_info=[l for l in currentSpe_STRINGcurrentSpePhyPPI_negPPI_info if  ((l[0],l[1]) not in currentSpe_string_score_dict) or (currentSpe_string_score_dict[(l[0],l[1])]==0)]
        currentSpe_STRINGcurrentSpePhyPPI_negPPI_info=[l for l in currentSpe_STRINGcurrentSpePhyPPI_negPPI_info if  ((l[0],l[1]) not in currentSpe_string_score_dict) ]

        print("len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_info)",len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_info))


        # In[ ]:





        # In[ ]:


        ## here we  limited frquence of protein in postive and negative dataset. 
        ## acturally only for negative now 
        ## futer post-filter step can be applied later , for example by Remove Ribosome protein or not 


        # final_testPP_info=get_PPIwithLimitedProFre(final_testPP_info)
        # print(len(final_testPP_info))



        currentSpe_STRINGcurrentSpePhyPPI_negPPI_info=get_PPIwithLimitedProFreByOr(currentSpe_STRINGcurrentSpePhyPPI_negPPI_info)
        print("len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_info):",len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_info))


        # In[ ]:


        currentSpe_STRINGcurrentSpePhyPPI_negPPI_ProFreDict=defaultdict(int)
        for p1,p2, _,_ in currentSpe_STRINGcurrentSpePhyPPI_negPPI_info:
            currentSpe_STRINGcurrentSpePhyPPI_negPPI_ProFreDict[p1] +=1
            currentSpe_STRINGcurrentSpePhyPPI_negPPI_ProFreDict[p2] +=1


        # In[ ]:


        print("len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_ProFreDict):",len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_ProFreDict))
        len([v for v in currentSpe_STRINGcurrentSpePhyPPI_negPPI_ProFreDict.values() if v >50])


        # In[ ]:


        plt.hist(currentSpe_STRINGcurrentSpePhyPPI_negPPI_ProFreDict.values(),bins=100)
        #plt.xscale('log')
        plt.yscale('log')
        plt.show()


        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:


        currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps=[(l[0],l[1]) for l in currentSpe_STRINGcurrentSpePhyPPI_posPPI_info]
        currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps=[(l[0],l[1]) for l in currentSpe_STRINGcurrentSpePhyPPI_negPPI_info]

        print("len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps):",len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps),len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps),len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps)+len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps))


        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # ## get paired MSA if pass certain creteria 

        # In[ ]:


        pairedMSA_unfiltered_folder


        # In[ ]:


        #%%time 

        pairedMSA_unfiltered_files =os.listdir(pairedMSA_unfiltered_folder) # here to check if paired MSA already been processed before 
        print("len(pairedMSA_unfiltered_files):",len(pairedMSA_unfiltered_files))
        print(pairedMSA_unfiltered_files[0:3])

        pairedMSA_unfiltered_pps=[f.split("and") for f in pairedMSA_unfiltered_files]
        pairedMSA_unfiltered_pps=[(p1,p2[0:-6]) for p1,p2 in pairedMSA_unfiltered_pps]

        pairedMSA_unfiltered_dict=dict([(f,1) for f in pairedMSA_unfiltered_pps])


        currentSpe_STRINGcurrentSpePhyPPI_posPPI_forPairMSA=[pp for pp in currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps if pp not in pairedMSA_unfiltered_dict]
        currentSpe_STRINGcurrentSpePhyPPI_posPPI_ArgForPairMSA=[(currentSpe_TaxID,currentSpe_hmmalign_path,p1,p2,Nf90_thres,pairedMSA_unfiltered_folder,pairedMSA_hhfilter_folder,pairedMSA_Nf90_folder) for p1,p2 in currentSpe_STRINGcurrentSpePhyPPI_posPPI_forPairMSA]

        currentSpe_STRINGcurrentSpePhyPPI_negPPI_forPairMSA=[pp for pp in currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps if pp not in pairedMSA_unfiltered_dict]
        currentSpe_STRINGcurrentSpePhyPPI_negPPI_ArgForPairMSA=[(currentSpe_TaxID,currentSpe_hmmalign_path,p1,p2,Nf90_thres,pairedMSA_unfiltered_folder,pairedMSA_hhfilter_folder,pairedMSA_Nf90_folder) for p1,p2 in currentSpe_STRINGcurrentSpePhyPPI_negPPI_forPairMSA]

        print("len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps),len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps)")
        print(len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps),len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps))
        print(len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_ArgForPairMSA),len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_ArgForPairMSA))

        # #run once #actually this is not necesseary  as function will check  if files is still existed 


        if len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_ArgForPairMSA)>0:
            #pool=mp.Pool(50)   #mp.Pool(160) 
            pool=get_context("spawn").Pool(50)
            get_pairedMSA_inOneRun_result=pool.map(get_pairedMSA_inOneRun,currentSpe_STRINGcurrentSpePhyPPI_posPPI_ArgForPairMSA)
            pool.close() 


            get_pairedMSA_inOneRun_result=[s for s in get_pairedMSA_inOneRun_result if s is not None]
            get_pairedMSA_inOneRun_result_frame=pd.DataFrame(get_pairedMSA_inOneRun_result,columns=["protein1","protein2","L1","L2","Nf90"])
            get_pairedMSA_inOneRun_result_frame.to_csv(pairedMSA_Nf90_csv,mode="a",
                                header=None,index=None,sep="\t")
            print(get_pairedMSA_inOneRun_result_frame.shape)



        if len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_ArgForPairMSA)>0:

            #pool=mp.Pool(50)   #mp.Pool(160) 
            pool=get_context("spawn").Pool(50)
            get_pairedMSA_inOneRun_result=pool.map(get_pairedMSA_inOneRun,currentSpe_STRINGcurrentSpePhyPPI_negPPI_ArgForPairMSA)
            pool.close() 


            get_pairedMSA_inOneRun_result=[s for s in get_pairedMSA_inOneRun_result if s is not None]
            get_pairedMSA_inOneRun_result_frame=pd.DataFrame(get_pairedMSA_inOneRun_result,columns=["protein1","protein2","L1","L2","Nf90"])
            get_pairedMSA_inOneRun_result_frame.to_csv(pairedMSA_Nf90_csv,mode="a",
                                header=None,index=None,sep="\t")
            print(get_pairedMSA_inOneRun_result_frame.shape)


        # CPU times: user 1.01 s, sys: 1.55 s, total: 2.56 s
        # Wall time: 5min 34s


        # In[ ]:


        #??? why so few ???


        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # In[ ]:





        # ## get postivte and negative ppi  second setp (after get all needed final paired MSA)
        # by  same protein ratio on two side of MSA 
        # 

        # In[ ]:


        # first need to update postive and negative ppi to only use pp that have large Nf90 value 
        pairedMSA_Nf90_frame = pd.read_csv(pairedMSA_Nf90_csv,header=None,index_col=None,sep="\t")
        print("pairedMSA_Nf90_frame.shape:",pairedMSA_Nf90_frame.shape)

        pairedMSA_Nf90_list=pairedMSA_Nf90_frame.iloc[:,[0,1]].values.tolist()
        pairedMSA_Nf90_dict=dict([((p1,p2),1) for p1, p2 in pairedMSA_Nf90_list])

        currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps=[pp for pp in currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps if pp in pairedMSA_Nf90_dict]
        currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps=[pp for pp in currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps if pp in pairedMSA_Nf90_dict]
        print("len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps):",len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps))
        print("len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps):",len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps))


        # In[ ]:


        #%%time 

        #heck paired MSA of homologus , if same proteins  from same species in both side

        if os.path.getsize(pairedMSA_sameProteinRatio_csv) == 0:
            pairedMSA_sameProteinRatio_dict=dict()
        else:
            pairedMSA_sameProteinRatio_frame = pd.read_csv(pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")
            pairedMSA_sameProteinRatio_list=pairedMSA_sameProteinRatio_frame.values.tolist()
            pairedMSA_sameProteinRatio_dict=dict([((p1,p2),r) for p1 , p2 , r in pairedMSA_sameProteinRatio_list])
            #print(pairedMSA_sameProteinRatio_frame.shape)

        currentSpe_STRINGcurrentSpePhyPPI_allPPI_pps=currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps+currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps
        currentSpe_STRINGcurrentSpePhyPPI_allPPI_ArgForSamePro=[(p1,p2,pairedMSA_Nf90_folder) for p1, p2, in currentSpe_STRINGcurrentSpePhyPPI_allPPI_pps if (p1, p2) not in pairedMSA_sameProteinRatio_dict]
        print("len(currentSpe_STRINGcurrentSpePhyPPI_allPPI_ArgForSamePro):",len(currentSpe_STRINGcurrentSpePhyPPI_allPPI_ArgForSamePro))

        if len(currentSpe_STRINGcurrentSpePhyPPI_allPPI_ArgForSamePro)>0:
            pool=mp.Pool(50)
            sameProtein_ratio_results=pool.map(getSameProteinRatio,currentSpe_STRINGcurrentSpePhyPPI_allPPI_ArgForSamePro)
            pool.close() 


            sameProtein_ratio_frame=pd.DataFrame(sameProtein_ratio_results,columns=["protein1","protein2","saemProtein_ratio"])
            sameProtein_ratio_frame.to_csv(pairedMSA_sameProteinRatio_csv,mode="a",
                                header=None,index=None,sep="\t")
            print(sameProtein_ratio_frame.shape)

        # CPU times: user 144 ms, sys: 639 ms, total: 783 ms
        # Wall time: 53.2 s


        # In[ ]:


        # remove protein pair with same prortein on two side of MSA 

        sameProtein_ratio_frame=pd.read_csv(pairedMSA_sameProteinRatio_csv,
                                header=None,index_col=None,sep="\t")
        print(sameProtein_ratio_frame.shape)
        large_sameProtein_ratio_frame=sameProtein_ratio_frame.loc[sameProtein_ratio_frame.iloc[:,2]>0,:]
        print(large_sameProtein_ratio_frame.shape)

        large_sameProtein_ratio_dict=dict([((p1,p2),r) for p1,p2,r in large_sameProtein_ratio_frame.values.tolist()])




        currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps=[pp for pp in currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps if pp not in large_sameProtein_ratio_dict]
        print("len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps):",len(currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps))


        currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps=[pp for pp in currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps if pp not in large_sameProtein_ratio_dict]
        print("len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps):",len(currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps))


        # In[ ]:




        currentSpe_STRINGcurrentSpePhyPPI_posPPI_beforeCoEvoComp_info=[(p1,p2, "P") for p1, p2 in currentSpe_STRINGcurrentSpePhyPPI_posPPI_pps]
        currentSpe_STRINGcurrentSpePhyPPI_negPPI_beforeCoEvoComp_info=[(p1,p2, "N") for p1, p2 in currentSpe_STRINGcurrentSpePhyPPI_negPPI_pps]

        currentSpe_STRINGcurrentSpePhyPPI_allPPI_beforeCoEvoComp_frame=pd.DataFrame(currentSpe_STRINGcurrentSpePhyPPI_posPPI_beforeCoEvoComp_info+currentSpe_STRINGcurrentSpePhyPPI_negPPI_beforeCoEvoComp_info,
                                                                                   columns=["STRING_ID1","STRING_ID2","benchmark_status"])
        currentSpe_STRINGcurrentSpePhyPPI_allPPI_beforeCoEvoComp_frame.to_csv(Benchmark_folder+"PPIInfoBeforeCoEvoComp.csv",
                                                                             header=True, index=False,sep="\t")






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
    except: 
        print("oops , some filter steps didnt pass, read log file; for  species and eggnog level:  ",currentSpe_TaxID ,";",EggNOG_maxLevel)
        message="oops , some filter steps didnt pass, read log file; for  species and eggnog level:  "+currentSpe_TaxID +";"+EggNOG_maxLevel
        with smtplib.SMTP(smtp_server, port) as server:
            server.ehlo()  # Can be omitted
            server.starttls(context=context)
            server.ehlo()  # Can be omitted
            server.login(sender_email, password)
            server.sendmail(sender_email, receiver_email, message)




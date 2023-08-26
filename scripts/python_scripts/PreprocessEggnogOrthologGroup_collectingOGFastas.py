#!/usr/bin/env python
# coding: utf-8

import argparse
import os 
import glob
import multiprocessing as mp
import pandas as pd
from collections import defaultdict


from create_singleMSA import combinedRBH_STRING1105


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
        # print(pro)
        pass
    else:
        return(RBH_HitsResults_Otholog_dic)

if __name__ == '__main__':
# run code in this way : python choose_orthologs_STRING11.05.py -o "fsd" -i "asdf" -m "asdfa"fsd asdf asdfa
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-fa','--currentSpe_fastaData', type=str, help='currentSpe_fastaData')
    parser.add_argument('-id','--currentSpe_TaxID', type=str, help='currentSpe_TaxID')
    parser.add_argument('-c','--currentSpe_currentMaxLevel_orthologs', type=str, help='currentSpe_currentMaxLevel_orthologs')
    parser.add_argument('-r','--redundant_proteins_csvFile', type=str, help='redundant_proteins_csvFile')
    parser.add_argument('-f','--newsingleMSA_RBH_OrthologousGroup_fileName', type=str, help='newsingleMSA_RBH_OrthologousGroup_fileName')
    parser.add_argument('-fa','--currentSpe_OrthologousGroup_Fa_path', type=str, help='currentSpe_OrthologousGroup_Fa_path')
    parser.add_argument('-log','--currentSpe_OrthologousGroup_Fa_logpath', type=str, help='currentSpe_OrthologousGroup_Fa_logpath') 
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    parser.add_argument('-ut','--code_utilities_folder', type=str, help='code_utilities_folder')

    args = parser.parse_args()
    currentSpe_fastaData=args.currentSpe_fastaData
    currentSpe_TaxID=args.currentSpe_TaxID
    currentSpe_currentMaxLevel_orthologs=args.currentSpe_currentMaxLevel_orthologs
    redundant_proteins_csvFile=args.redundant_proteins_csvFile
    newsingleMSA_RBH_OrthologousGroup_fileName=args.newsingleMSA_RBH_OrthologousGroup_fileName
    origSTRINGBacteriaProSeqPath=args.origSTRINGBacteriaProSeqPath
    currentSpe_OrthologousGroup_Fa_path=args.currentSpe_OrthologousGroup_Fa_path
    currentSpe_OrthologousGroup_Fa_logpath=args.currentSpe_OrthologousGroup_Fa_logpath
    mp_task_nums=int(args.mp_task_nums)
    code_utilities_folder=args.code_utilities_folder
    
    
    origSTRINGBacteriaFiles=glob.glob(os.path.join(origSTRINGBacteriaProSeqPath,"*.fa"))
    print("len(origSTRINGBacteriaFiles):",len(origSTRINGBacteriaFiles))
    
    # get new othologous group
    print("args.currentSpe_currentMaxLevel_orthologs",args.currentSpe_currentMaxLevel_orthologs)
    generate_ogs_files=glob.glob(os.path.join(currentSpe_currentMaxLevel_orthologs,"*.orthologs")) # here currentSpe_currentMaxLevel_orthologs from nextflow output channel witl remove "/" by defaults, thus use here. 
    print("len(generate_ogs_files):",len(generate_ogs_files))


    # check how many of them orthologs files are not empty  or morn then a threshold 
    # here filter by number of orthologs in single MSA 
    generate_ogs_ArgForCombinedRBH=[(f,len(origSTRINGBacteriaFiles) )for f in generate_ogs_files]
    pool = mp.Pool(mp_task_nums) #better set this as a nextflow param as well ?
    combinedRBH_results = pool.map(combinedRBH_STRING1105, generate_ogs_ArgForCombinedRBH)
    pool.close()

    print("len(combinedRBH_results):",len(combinedRBH_results))
    combinedRBH_results=[c for c in combinedRBH_results if c is not None]
    print("len(combinedRBH_results):",len(combinedRBH_results))
    combinedRBH_results_flatten=pd.concat(combinedRBH_results)
    print("combinedRBH_results_flatten:",len(combinedRBH_results_flatten))


    #remove reduant proteins by checking seq similarities 
    currentSpe_redundant_pros=pd.read_csv(redundant_proteins_csvFile,
                                     header=0,index_col=None)
    print("currentSpe_redundant_pros.shape:",currentSpe_redundant_pros.shape)
    currentSpe_redundant_pros_list=list(currentSpe_redundant_pros['pro_name'])

    NonRedundantIdx=[False if pro in currentSpe_redundant_pros_list else True for pro in combinedRBH_results_flatten.iloc[:,0].tolist()]


    combinedRBH_results_frame=combinedRBH_results_flatten.rename(columns={0:"query_protein",1:"RBH_hit",})
    print("combinedRBH_results_frame.shape:",combinedRBH_results_frame.shape)
    combinedRBH_results_frame_NonRedundantIdx=combinedRBH_results_frame.loc[NonRedundantIdx,:]
    print("combinedRBH_results_frame_NonRedundantIdx.shape:",combinedRBH_results_frame_NonRedundantIdx.shape)
    combinedRBH_results_frame_NonRedundantIdx.head()

    Unique_currentSpe_Pros=list(set(combinedRBH_results_frame_NonRedundantIdx['query_protein']))
    #Unique_Ecoli_Pros=list(set(Ecoli_within.iloc[:,0]))

    print("len(Unique_currentSpe_Pros),Unique_currentSpe_Pros[0:3]",len(Unique_currentSpe_Pros),Unique_currentSpe_Pros[0:3]) #check this value



    # this is original code, in nextflow , dont need to check file existence since they are in current process working directory , so adapt it 
    # newsingleMSA_RBH_OrthologousGroup_fileName=currentSpeMiddleDataPath+"newsingleMSA_RBH_OrthologousGroup.csv"
    # if not os.path.exists(newsingleMSA_RBH_OrthologousGroup_fileName):
    #     pool = mp.Pool(30) # 
    #     results = pool.map(getCurrentSpeRBHdict, Unique_currentSpe_Pros) # here many damian generalted files are emply 
    #     pool.close() 

    #     results=[temp for temp in results if temp is not None]
    #     print("len(results):",len(results))


    #     # only run once as its append mode 
    #     for RBH_HitsResults_Otholog_dic_temp in results:
    #         with open(newsingleMSA_RBH_OrthologousGroup_fileName,"a") as file:
    #             for key in RBH_HitsResults_Otholog_dic_temp.keys():
    #                 file.write("%s,%s\n" % (key, ",".join(RBH_HitsResults_Otholog_dic_temp[key]))) 



    pool = mp.Pool(mp_task_nums) # 
    results = pool.map(getCurrentSpeRBHdict, Unique_currentSpe_Pros) # here many damian generalted files are emply 
    pool.close() 

    results=[temp for temp in results if temp is not None]
    print("len(results):",len(results))

    # only run once as its append mode , nextflow model, eachtime , if process restart, the new files is created , so no problem, compare file size of multiple run 
    for RBH_HitsResults_Otholog_dic_temp in results:
        with open(newsingleMSA_RBH_OrthologousGroup_fileName,"a") as file: # normally this command could just create a new file, but in nextflow, this command somehow failed , why ?,oh the problem is that i didnt create the folder 
            for key in RBH_HitsResults_Otholog_dic_temp.keys():
                file.write("%s,%s\n" % (key, ",".join(RBH_HitsResults_Otholog_dic_temp[key]))) 


                
                

    #notice following code could be simply all implemented in bash script 
    #here for we use python to wrap bash code because Im familiar with parallel processing in python 
    newsingleMSA_RBH_OrthologousGroup_frame=pd.read_csv(newsingleMSA_RBH_OrthologousGroup_fileName,
                                                        header=None,index_col=None, sep=" ")
    print(newsingleMSA_RBH_OrthologousGroup_frame.shape)
    newsingleMSA_RBH_OrthologousGroup_OGidx=[idx+1 for idx in range(newsingleMSA_RBH_OrthologousGroup_frame.shape[0])]

    newsingleMSA_RBH_OrthologousGroup_OGidx_ArgForFaidx=[(code_utilities_folder,currentSpe_TaxID,OG_idx,currentSpe_fastaData,origSTRINGBacteriaProSeqPath,currentSpe_OrthologousGroup_Fa_path,currentSpe_OrthologousGroup_Fa_logpath,newsingleMSA_RBH_OrthologousGroup_fileName) for OG_idx in newsingleMSA_RBH_OrthologousGroup_OGidx]
    pool=mp.Pool(mp_task_nums)  # here when we set it as 50, there are more i.o wait time from top command ; 
    pool.map(fun_newSingleMSA_EggNOG_OrthologousGroup_faidx,newsingleMSA_RBH_OrthologousGroup_OGidx_ArgForFaidx)
    pool.close() 

        
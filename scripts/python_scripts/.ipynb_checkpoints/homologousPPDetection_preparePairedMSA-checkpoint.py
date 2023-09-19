
from __future__ import division
import argparse
import re
import multiprocessing as mp
from multiprocessing import get_context
import glob
import pandas as pd
import sys
import os
import pickle 


from  create_pairedMSA import get_pairedMSA_inOneRun
from create_pairedMSA import  getSameProteinRatio


if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument('-q','--Query_tuple', type=str, help='Query_tuple')
    parser.add_argument('-s','--Subject_tupleList', type=str, help='Subject_tupleList')
    parser.add_argument('-f','--newSTRING_rootFolder', type=str, help='newSTRING_rootFolder')
    parser.add_argument('-c','--CoEvo_data_folder', type=str, help='CoEvo_data_folder')
    parser.add_argument('-mb','--homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path', type=str, help='homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path')

    parser.add_argument('-nf90','--Nf90_thres', type=str, help='Nf90_thres')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    
    args = parser.parse_args()
    Query_tuple=args.Query_tuple
    Subject_tupleList=args.Subject_tupleList
    newSTRING_rootFolder=args.newSTRING_rootFolder
    CoEvo_data_folder=args.CoEvo_data_folder
    homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path=args.homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path
    
    Nf90_thres=int(args.Nf90_thres)
    mp_task_nums=int(args.mp_task_nums)
    
    Query_tuple=Query_tuple.split("_")
    Query_speID=Query_tuple[1]
    Subject_tupleList=Subject_tupleList.split("_")
    Subject_tupleList=[(Subject_tupleList[i],Subject_tupleList[i+1]) for i in range(0,len(Subject_tupleList)-1,2)]
    
    
    
    query_EggNOG_maxLevel,query_TaxID=Query_tuple
    
    with open(homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_BestHomologous_listDict.pickle", 'rb') as handle:
        Query2Subject_QueSpeAllPPI_BestHomologous_listDict=pickle.load(handle)  



    for currentSubject_EggNOG_maxLevel,currentSubject_TaxID in Subject_tupleList:
        print(currentSubject_EggNOG_maxLevel,currentSubject_TaxID)
        # adjust these three folder to fit to nextlfow convention
        # newSTRING_rootFolder="/mnt/mnemo6/tao/PPI_Coevolution/STRING_data_11.5/"
        # CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING11.5/"
        #     input_root_folder=CoEvo_data_folder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_eggNOGfilteredData/"

        #here is che coevo_data_folder is the nextflow  temperaoly folder , 
        #then the created folder is also in the nextflow temperaroly folder ?
        
        #if without prefix, its for subject species 
        input_root_folder=CoEvo_data_folder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_eggNOGfilteredData/"
        pairedMSA_unfiltered_folder=input_root_folder+"pair_MSA_unfiltered_PasteAlign/"
        pairedMSA_hhfilter_folder=input_root_folder+"pair_MSA_hhfilter_PasteAlign/"
        pairedMSA_Nf90_folder = input_root_folder+"pair_MSA_Nf90_PasteAlign/"

        # DCA_coevolutoin_path=input_root_folder+"coevolutoin_result_DCA/"
        # MI_coevolutoin_path=input_root_folder+"coevolutoin_result_MI/"

        #this two folder denition is inherated from previous pipeline
        #previous  we run the the created pairedMSA pipeline for subject species independently as we did for query speceis
        #thats why we have these two folders. but since now i first created paired MSA for all query pp
        # then we created for all homoglos subject pp of all query pps, we dont need this two step. but no harm to leave them here
        original_pairedMSA_Nf90_csv = input_root_folder+"pair_MSA_Nf90.csv"
        original_pairedMSA_sameProteinRatio_csv = input_root_folder+"sameProteinRatio.csv"

        Query_prefix="BestHomologousPPFor"+Query_tuple[1]+"AtEggNOGmaxLevel"+Query_tuple[0]+"_"


        pairedMSA_Nf90_csv = input_root_folder+Query_prefix+"pair_MSA_Nf90.csv"
        pairedMSA_sameProteinRatio_csv = input_root_folder+Query_prefix+"sameProteinRatio.csv"

        Benchmark_folder=input_root_folder+Query_prefix+"AllPPI_Benchmark/" #

        currentSubject_hmmalign_path=newSTRING_rootFolder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_newSingleMSA_hmmalign_removeGaps_keepGapPos/"
        currentSpeMSAGapsFilteringMetaFolder=newSTRING_rootFolder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_MSAGapsFiltering/"


        if not os.path.exists(Benchmark_folder):
            os.makedirs(Benchmark_folder)

#         if not os.path.exists(DCA_coevolutoin_path):
#             os.makedirs(DCA_coevolutoin_path)

#         if not os.path.exists(MI_coevolutoin_path):
#             os.makedirs(MI_coevolutoin_path)

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

        if not os.path.exists(original_pairedMSA_Nf90_csv):
            with open(original_pairedMSA_Nf90_csv, 'w') as fp:
                pass

        if not os.path.exists(original_pairedMSA_sameProteinRatio_csv):
            with open(original_pairedMSA_sameProteinRatio_csv, 'w') as fp:
                pass
       


        currentSubject_Query2Subject_QueSpeAllPPI_BestHomologous_dict=Query2Subject_QueSpeAllPPI_BestHomologous_listDict[(currentSubject_EggNOG_maxLevel,currentSubject_TaxID)]

        print("len(currentQuery_Subject2Query_QueSpeAllPPI_BestHomologous_dict):",len(currentSubject_Query2Subject_QueSpeAllPPI_BestHomologous_dict))
        print(list(currentSubject_Query2Subject_QueSpeAllPPI_BestHomologous_dict.items())[0:3])


        Current_Subject_homologousPPs=list()
        for que_pp , best_subject_pp in currentSubject_Query2Subject_QueSpeAllPPI_BestHomologous_dict.items():
            Current_Subject_homologousPPs.append(best_subject_pp)
        print("len(Current_Subject_homologousPPs):",len(Current_Subject_homologousPPs))
        Current_Subject_homologousPPs=list(set(Current_Subject_homologousPPs))
        print("len(set(Current_Subject_homologousPPs)):",len(Current_Subject_homologousPPs))



        with open(currentSpeMSAGapsFilteringMetaFolder+'fasta_protein_lens_dict.pickle', 'rb') as handle:
            fasta_protein_lens=pickle.load(handle)

        with open(currentSpeMSAGapsFilteringMetaFolder+'fasta_protein_Nf90s_dict.pickle', 'rb') as handle:
            fasta_protein_Nf90s=pickle.load(handle)




        currentSubject_final_MSAs=glob.glob(currentSubject_hmmalign_path+"*")
        print("len(currentSubject_final_MSAs:",len(currentSubject_final_MSAs))
        print(currentSubject_final_MSAs[0])

        currentSubject_final_pros=[os.path.basename(f) for f in currentSubject_final_MSAs]
        currentSubject_final_pros=[p[:-3] for p in currentSubject_final_pros]
        print(currentSubject_final_pros[0])
        currentSubject_final_pro_dict=dict([(p,1)for p in currentSubject_final_pros])



        pairedMSA_unfiltered_dict=dict()

        currentSubject_PPI_forPairMSA=[pp for pp in Current_Subject_homologousPPs if pp not in pairedMSA_unfiltered_dict]
        currentSubject_PPI_ArgForPairMSA=[(currentSubject_TaxID,currentSubject_hmmalign_path,p1,p2,Nf90_thres,pairedMSA_unfiltered_folder,pairedMSA_hhfilter_folder,pairedMSA_Nf90_folder) for p1,p2 in currentSubject_PPI_forPairMSA]


        print("len(Current_Subject_homologousPPs):",len(Current_Subject_homologousPPs))
        print("len(currentSubject_PPI_ArgForPairMSA):",len(currentSubject_PPI_ArgForPairMSA))




        if len(currentSubject_PPI_ArgForPairMSA)>0:

            pool=mp.Pool(mp_task_nums)   #mp.Pool(160) 
            get_pairedMSA_inOneRun_result=pool.map(get_pairedMSA_inOneRun,currentSubject_PPI_ArgForPairMSA)
            pool.close() 


            get_pairedMSA_inOneRun_result=[s for s in get_pairedMSA_inOneRun_result if s is not None]
            get_pairedMSA_inOneRun_result_frame=pd.DataFrame(get_pairedMSA_inOneRun_result,columns=["protein1","protein2","L1","L2","Nf90"])
            get_pairedMSA_inOneRun_result_frame.to_csv(pairedMSA_Nf90_csv,mode="a",
                                header=None,index=None,sep="\t")
            print("get_pairedMSA_inOneRun_result_frame.shape:",get_pairedMSA_inOneRun_result_frame.shape)




        # ## get postivte and negative ppi  second setp (after get all needed final paired MSA)
        # by  same protein ratio on two side of MSA 

        # first need to update postive and negative ppi to only use pp that have large Nf90 value 
        if os.path.getsize(original_pairedMSA_Nf90_csv) == 0:
            original_pairedMSA_Nf90_frame=pd.DataFrame()
        else:
            original_pairedMSA_Nf90_frame = pd.read_csv(original_pairedMSA_Nf90_csv,header=None,index_col=None,sep="\t")
            print("original_pairedMSA_Nf90_frame.shape:",original_pairedMSA_Nf90_frame.shape)
        
        pairedMSA_Nf90_frame = pd.read_csv(pairedMSA_Nf90_csv,header=None,index_col=None,sep="\t")
        print("pairedMSA_Nf90_frame.shape:",pairedMSA_Nf90_frame.shape)
        
        pairedMSA_Nf90_frame=pd.concat([original_pairedMSA_Nf90_frame,pairedMSA_Nf90_frame,]) 
        print("pairedMSA_Nf90_frame.shape:",pairedMSA_Nf90_frame.shape)
        pairedMSA_Nf90_list=pairedMSA_Nf90_frame.iloc[:,[0,1]].values.tolist()
        pairedMSA_Nf90_dict=dict([((p1,p2),1) for p1, p2 in pairedMSA_Nf90_list])

        print(len(pairedMSA_Nf90_dict))
        print(len(Current_Subject_homologousPPs))
        Current_Subject_homologousPPs=[pp for pp in Current_Subject_homologousPPs if pp in pairedMSA_Nf90_dict]
        print("len(Current_Subject_homologousPPs):",len(Current_Subject_homologousPPs))





        #heck paired MSA of homologus , if same proteins  from same species in both side

        if os.path.getsize(original_pairedMSA_sameProteinRatio_csv) == 0:
            original_pairedMSA_sameProteinRatio_frame=pd.DataFrame()
        else:
            original_pairedMSA_sameProteinRatio_frame = pd.read_csv(original_pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")
            
        if os.path.getsize(pairedMSA_sameProteinRatio_csv) == 0:
            pairedMSA_sameProteinRatio_frame=pd.DataFrame()
        else:
            pairedMSA_sameProteinRatio_frame = pd.read_csv(pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")
            
            

        pairedMSA_sameProteinRatio_frame=pd.concat([original_pairedMSA_sameProteinRatio_frame,pairedMSA_sameProteinRatio_frame])
        pairedMSA_sameProteinRatio_list=pairedMSA_sameProteinRatio_frame.values.tolist()
        pairedMSA_sameProteinRatio_dict=dict([((p1,p2),r) for p1 , p2 , r in pairedMSA_sameProteinRatio_list])
        #print(pairedMSA_sameProteinRatio_frame.shape)


        Current_Subject_homologousPPs_ArgForSamePro=[(p1,p2,pairedMSA_Nf90_folder) for p1, p2, in Current_Subject_homologousPPs if (p1, p2) not in pairedMSA_sameProteinRatio_dict]
        print("len(Current_Subject_homologousPPs_ArgForSamePro):",len(Current_Subject_homologousPPs_ArgForSamePro))


        # In[ ]:


        if len(Current_Subject_homologousPPs_ArgForSamePro)>0:
            pool=mp.Pool(mp_task_nums)
            sameProtein_ratio_results=pool.map(getSameProteinRatio,Current_Subject_homologousPPs_ArgForSamePro)
            pool.close() 


            sameProtein_ratio_frame=pd.DataFrame(sameProtein_ratio_results,columns=["protein1","protein2","saemProtein_ratio"])
            sameProtein_ratio_frame.to_csv(pairedMSA_sameProteinRatio_csv,mode="a",
                                header=None,index=None,sep="\t")
            print("sameProtein_ratio_frame.shape:",sameProtein_ratio_frame.shape)

        # CPU times: user 144 ms, sys: 639 ms, total: 783 ms
        # Wall time: 53.2 s



        # remove protein pair with same prortein on two side of MSA 

        sameProtein_ratio_frame = pd.read_csv(pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")

        if os.path.getsize(original_pairedMSA_sameProteinRatio_csv) == 0:
            original_sameProtein_ratio_frame=pd.DataFrame()
        else:
            original_sameProtein_ratio_frame = pd.read_csv(original_pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")


        sameProtein_ratio_frame=pd.concat([original_sameProtein_ratio_frame,sameProtein_ratio_frame])

        print("sameProtein_ratio_frame.shape:",sameProtein_ratio_frame.shape)

        large_sameProtein_ratio_frame=sameProtein_ratio_frame.loc[sameProtein_ratio_frame.iloc[:,2]>0,:]
        print("large_sameProtein_ratio_frame.shape:",large_sameProtein_ratio_frame.shape)

        large_sameProtein_ratio_dict=dict([((p1,p2),r) for p1,p2,r in large_sameProtein_ratio_frame.values.tolist()])



        Current_Subject_homologousPPs=[pp for pp in Current_Subject_homologousPPs if pp not in large_sameProtein_ratio_dict]
        print(len(Current_Subject_homologousPPs))



        Current_Subject_homologousPPs_beforeDCAcomputation_file=Benchmark_folder+"Current_Subject_BestHomologousPPs_beforeDCAcomputation.pickle"
        if not os.path.exists(Current_Subject_homologousPPs_beforeDCAcomputation_file):
            with open(Current_Subject_homologousPPs_beforeDCAcomputation_file, 'wb') as handle:
                    pickle.dump(Current_Subject_homologousPPs,handle)



        
        


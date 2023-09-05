
from __future__ import division
import argparse
import re
import multiprocessing as mp
from multiprocessing import get_context
import glob
import pandas as pd
import sys
import os






if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument('-q','--query_EggNOG_maxLevel', type=str, help='query_EggNOG_maxLevel')
    parser.add_argument('-d','--query_TaxID', type=str, help='query_TaxID')

    parser.add_argument('-s','--currentSubect_EggNOG_maxLevel', type=str, help='currentSubject_EggNOG_maxLevel')
    parser.add_argument('-i','--currentSubject_TaxID', type=str, help='currentSubject_TaxID')



    args = parser.parse_args()
    query_EggNOG_maxLevel=args.query_EggNOG_maxLevel
    query_TaxID=args.query_TaxID

    currentSubject_EggNOG_maxLevel=args.currentSubect_EggNOG_maxLevel
    currentSubject_TaxID=args.currentSubject_TaxID


    CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING11.5/"


    Query_tuple=(query_EggNOG_maxLevel, query_TaxID)





    Query_input_root_folder=CoEvo_data_folder+Query_tuple[1]+"_EggNOGmaxLevel"+Query_tuple[0]+"_eggNOGfilteredData/"
    Query_Benchmark_folder=Query_input_root_folder+"STRINGPhyBalancePhyla_Benchmark/"



    
    #if without prefix, its for subject species 
    input_root_folder=CoEvo_data_folder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_eggNOGfilteredData/"
    pairedMSA_unfiltered_folder=input_root_folder+"pair_MSA_unfiltered_PasteAlign/"
    pairedMSA_hhfilter_folder=input_root_folder+"pair_MSA_hhfilter_PasteAlign/"
    pairedMSA_Nf90_folder = input_root_folder+"pair_MSA_Nf90_PasteAlign/"

    DCA_coevolutoin_path=input_root_folder+"coevolutoin_result_DCA/"
    MI_coevolutoin_path=input_root_folder+"coevolutoin_result_MI/"

    original_pairedMSA_Nf90_csv = input_root_folder+"pair_MSA_Nf90.csv"
    original_pairedMSA_sameProteinRatio_csv = input_root_folder+"sameProteinRatio.csv"

    Query_prefix="BestHomologousPPFor"+Query_tuple[1]+"AtEggNOGmaxLevel"+Query_tuple[0]+"_"
    pairedMSA_Nf90_csv = input_root_folder+Query_prefix+"pair_MSA_Nf90.csv"
    pairedMSA_sameProteinRatio_csv = input_root_folder+Query_prefix+"sameProteinRatio.csv"


    changtge this to a file a later 
    Benchmark_folder=input_root_folder+Query_prefix+"STRINGPhyBalancePhyla_Benchmark/" #


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



    newSTRING_rootFolder="/mnt/mnemo6/tao/PPI_Coevolution/STRING_data_11.5/"
    currentSubject_hmmalign_path=newSTRING_rootFolder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_newSingleMSA_hmmalign_removeGaps_keepGapPos/"

    currentSubjectMiddleDataPath=newSTRING_rootFolder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_MiddleData/"







    with open("/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING11.5/511145_EggNOGmaxLevel1224_eggNOGfilteredData/STRINPhyPPI_Benchmark/NameUnsorted_Subject2Query_SubSpeAllPPI_BestHomologous_ignoreQueryDCA_dict_listDict.pickle", 'rb') as handle:
        Query2Subject_SubSpePPHighDCA_BestHomologous_ignoreQueryDCA_dict_listDict=pickle.load(handle)



    currentSubject_Query2Subject_SubSpePPHighDCA_BestHomologous_ignoreQueryDCA_dict=Query2Subject_SubSpePPHighDCA_BestHomologous_ignoreQueryDCA_dict_listDict[(currentSubject_EggNOG_maxLevel,currentSubject_TaxID)]



    print("len(currentQuery_Subject2Query_SubSpePPHighDCA_BestHomologous_ignoreQueryDCA_dict):",len(currentSubject_Query2Subject_SubSpePPHighDCA_BestHomologous_ignoreQueryDCA_dict))
    print(list(currentSubject_QuerySubject_SubSpePPHighDCA_BestHomologous_ignoreQueryDCA_dict.items())[0:3])



    # In[ ]:


    count=0
    Current_Subject_homologousPPs=list()
    #for sub_pp , query_pps in currentQuery_Subject2Query_SubSpePPHighDCA_BestHomologous_ignoreQueryDCA_dict.items():
    for que_pp , best_subject_pp in currentSubject_Query2Subject_SubSpePPHighDCA_BestHomologous_ignoreQueryDCA_dict.items():
        #sub_pp_info=HighDCA_Namesorted_Subject_allPPI_allInfo_dict[sub_pp]
        #sub_pp_status=sub_pp_info[0]

        #if sub_pp_status=="N":
        count +=1
        #Current_Query_homologousPPs.append(query_pps)
        Current_Query_homologousPPs.append(best_subject_pp)
        #print(sub_pp_info)

    print(count)



    # In[ ]:


    print("len(Current_Subject_homologousPPs),len(set(Current_Subject_homologousPPs)):",len(Current_Subject_homologousPPs),len(set(Current_Subject_homologousPPs)))


    # In[ ]:


    Current_Subject_homologousPPs=list(set(Current_Subject_homologousPPs))
    print(len(Current_Subject_homologousPPs))






    # In[ ]:


    Nf90_thres=16

    with open(currentSubjectMiddleDataPath+'fasta_protein_lens_dict.pickle', 'rb') as handle:
        fasta_protein_lens=pickle.load(handle)

    with open(currentSubjectMiddleDataPath+'fasta_protein_Nf90s_dict.pickle', 'rb') as handle:
        fasta_protein_Nf90s=pickle.load(handle)




    currentSubject_final_MSAs=glob.glob(currentSubject_hmmalign_path+"*")
    print("len(currentSubject_final_MSAs:",len(currentSubject_final_MSAs))
    print(currentSubject_final_MSAs[0])

    currentSubject_final_pros=[os.path.basename(f) for f in currentSubject_final_MSAs]
    currentSubject_final_pros=[p[:-3] for p in currentSubject_final_pros]
    print(currentSubject_final_pros[0])

    currentSubject_final_pro_dict=dict([(p,1)for p in currentSubject_final_pros])



    pairedMSA_unfiltered_dict=dict()

    currentSubject_STRINGPhyBalancePhylaPPI_forPairMSA=[pp for pp in Current_Subject_homologousPPs if pp not in pairedMSA_unfiltered_dict]
    currentSubject_STRINGPhyBalancePhylaPPI_ArgForPairMSA=[(currentSubject_TaxID,currentSubject_hmmalign_path,p1,p2,Nf90_thres,pairedMSA_unfiltered_folder,pairedMSA_hhfilter_folder,pairedMSA_Nf90_folder) for p1,p2 in currentSubject_STRINGPhyBalancePhylaPPI_forPairMSA]


    print("len(Current_Subject_homologousPPs):",len(Current_Subject_homologousPPs))
    print("len(currentSubject_STRINGPhyBalancePhylaPPI_ArgForPairMSA):",len(currentSubject_STRINGPhyBalancePhylaPPI_ArgForPairMSA))




    if len(currentSubject_STRINGPhyBalancePhylaPPI_ArgForPairMSA)>0:

        now = datetime.now()
        print("paired MSA one run starts now =", now)

        pool=mp.Pool(30)   #mp.Pool(160) 
        get_pairedMSA_inOneRun_result=pool.map(get_pairedMSA_inOneRun,currentSubject_STRINGPhyBalancePhylaPPI_ArgForPairMSA)
        pool.close() 


        get_pairedMSA_inOneRun_result=[s for s in get_pairedMSA_inOneRun_result if s is not None]
        get_pairedMSA_inOneRun_result_frame=pd.DataFrame(get_pairedMSA_inOneRun_result,columns=["protein1","protein2","L1","L2","Nf90"])
        get_pairedMSA_inOneRun_result_frame.to_csv(pairedMSA_Nf90_csv,mode="a",
                            header=None,index=None,sep="\t")
        print("get_pairedMSA_inOneRun_result_frame.shape:",get_pairedMSA_inOneRun_result_frame.shape)

    now = datetime.now()
    print("paired MSA one run ends now =", now)

  

    # ## get postivte and negative ppi  second setp (after get all needed final paired MSA)
    # by  same protein ratio on two side of MSA 

    # first need to update postive and negative ppi to only use pp that have large Nf90 value 
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

    if os.path.getsize(pairedMSA_sameProteinRatio_csv) == 0:
        pairedMSA_sameProteinRatio_dict=dict()
    else:
        pairedMSA_sameProteinRatio_frame = pd.read_csv(pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")
        original_pairedMSA_sameProteinRatio_frame = pd.read_csv(original_pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")
        pairedMSA_sameProteinRatio_frame=pd.concat([original_pairedMSA_sameProteinRatio_frame,pairedMSA_sameProteinRatio_frame])
        pairedMSA_sameProteinRatio_list=pairedMSA_sameProteinRatio_frame.values.tolist()
        pairedMSA_sameProteinRatio_dict=dict([((p1,p2),r) for p1 , p2 , r in pairedMSA_sameProteinRatio_list])
        #print(pairedMSA_sameProteinRatio_frame.shape)


    Current_Subject_homologousPPs_ArgForSamePro=[(p1,p2,pairedMSA_Nf90_folder) for p1, p2, in Current_Subject_homologousPPs if (p1, p2) not in pairedMSA_sameProteinRatio_dict]
    print("len(Current_Subject_homologousPPs_ArgForSamePro):",len(Current_Subject_homologousPPs_ArgForSamePro))


    # In[ ]:


    if len(Current_Subject_homologousPPs_ArgForSamePro)>0:
        pool=mp.Pool(50)
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

    original_sameProtein_ratio_frame = pd.read_csv(original_pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")


    sameProtein_ratio_frame=pd.concat([original_sameProtein_ratio_frame,sameProtein_ratio_frame])

    print("sameProtein_ratio_frame.shape:",sameProtein_ratio_frame.shape)

    large_sameProtein_ratio_frame=sameProtein_ratio_frame.loc[sameProtein_ratio_frame.iloc[:,2]>0,:]
    print("large_sameProtein_ratio_frame.shape:",large_sameProtein_ratio_frame.shape)

    large_sameProtein_ratio_dict=dict([((p1,p2),r) for p1,p2,r in large_sameProtein_ratio_frame.values.tolist()])



    Current_Query_homologousPPs=[pp for pp in Current_Query_homologousPPs if pp not in large_sameProtein_ratio_dict]
    print(len(Current_Query_homologousPPs))



    Current_Subject_homologousPPs_beforeDCAcomputation_file=Benchmark_folder+"Current_Subject_BestHomologousPPs_beforeDCAcomputation.pickle"
    if not os.path.exists(Current_Subject_homologousPPs_beforeDCAcomputation_file):
        with open(Current_Subject_homologousPPs_beforeDCAcomputation_file, 'wb') as handle:
                pickle.dump(Current_Subject_homologousPPs,handle)



        
        


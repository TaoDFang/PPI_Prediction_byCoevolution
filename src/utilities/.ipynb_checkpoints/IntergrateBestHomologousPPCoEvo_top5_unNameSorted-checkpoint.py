import pandas as pd
import numpy as np
import copy
import os
import sys
import pickle
import multiprocessing as mp 


from collections import defaultdict
import matplotlib.pyplot as plt 


from collect_topCoEvos import get_topRanking_CoEvo_file
from IntergrateBestHomologousPPCoEvo_unNameSorted import getMetaFrame_FullBestHomologousPP



def get_SubjectInfo_top5DCAs_FullBestHomologousPP(Query_tuple,Subject_tupleList,
                                                STRING_version="11.5",
                                                coevo_suffix="_pydcaFNAPC_array",
                                                benchmark_suffix="STRINPhyPPI_Benchmark/",
                                            BestHomologousPP_filePrefix="BestHomologousPP",
                                            CoEvo_data_folder="CoEvo_data_STRING11.5/",
                                                topDCA_num=5,
                                                returnDic=False,
                                                overwrite=True,
                                                use_multiprocessing=True,
                                               ):
    # for quey data 
    Subect_BestHomologousPP_prefix="BestHomologousPPFor"+Query_tuple[1]+"AtEggNOGmaxLevel"+Query_tuple[0]+"_"

    
    topDCA_idx=[3*i for i in range(topDCA_num)]
    
    Subject_allPPI_allInfo_dict_listDict=dict()

    for Subject_phylum, Subject_speID in Subject_tupleList:
        print(Subject_phylum,Subject_speID)

        Subject_allPPI_allInfo_frame=getMetaFrame_FullBestHomologousPP(Query_tuple=Query_tuple,EggNOG_maxLevel=Subject_phylum,currentSpe_TaxID=Subject_speID,
                                                              STRING_Version=STRING_version,benchmark_suffix=benchmark_suffix,smallPhylum="",BestHomologousPP_filePrefix=BestHomologousPP_filePrefix,CoEvo_data_folder=CoEvo_data_folder)
        Subject_allPPI_info=Subject_allPPI_allInfo_frame.loc[:,["STRING_ID1","STRING_ID2","len1","len2"]].values.tolist()
        
        
        smallPhylum=""
        input_root_folder=CoEvo_data_folder+smallPhylum+Subject_speID+"_EggNOGmaxLevel"+Subject_phylum+"_eggNOGfilteredData/"
        BestHomologousPP_Benchmark_folder=input_root_folder+Subect_BestHomologousPP_prefix+benchmark_suffix
        print(f"BestHomologousPP_Benchmark_folder:{BestHomologousPP_Benchmark_folder}")
        
        if coevo_suffix=="_pydcaFNAPC_array":
            DCA_coevolutoin_path=input_root_folder+"coevolutoin_result_DCA/"
            topRanking_pydcaFNAPC_file=BestHomologousPP_Benchmark_folder+BestHomologousPP_filePrefix+"_topRanking_pydcaFNAPC_frame.csv" 
        
            
        if coevo_suffix=="_alphafold_prob12":
            colab_output_path="/mnt/mnemo6/tao/colabfold_PPI_Coevolution/CoEvo_data_STRING"+STRING_version+"/"
            DCA_coevolutoin_path=os.path.join(colab_output_path,f"EggNog_{Subject_speID}_EggNOGmaxLeve{Subject_phylum}/")
            topRanking_pydcaFNAPC_file=BestHomologousPP_Benchmark_folder+BestHomologousPP_filePrefix+"_topRanking_alphafoldProb12_frame.csv" 
            
        if coevo_suffix=="_alphafold_inverseMinDist":
            colab_output_path="/mnt/mnemo6/tao/colabfold_PPI_Coevolution/CoEvo_data_STRING"+STRING_version+"/"
            DCA_coevolutoin_path=os.path.join(colab_output_path,f"EggNog_{Subject_speID}_EggNOGmaxLeve{Subject_phylum}/")
            topRanking_pydcaFNAPC_file=BestHomologousPP_Benchmark_folder+BestHomologousPP_filePrefix+"_topRanking_alphafoldinverseMinDist_frame.csv" 
        
        
        print(f"DCA_coevolutoin_path:{DCA_coevolutoin_path}")
        print(f"topRanking_pydcaFNAPC_file:{topRanking_pydcaFNAPC_file}")
        
        top_pydcaFNAPC_dict=get_topRanking_CoEvo_file(topRanking_CoEvo_file=topRanking_pydcaFNAPC_file,
                                                       coevolutoin_path=DCA_coevolutoin_path,
                                                       coevo_suffix=coevo_suffix,
                                                       allPPI_info=Subject_allPPI_info, 
                                                       returnDic=returnDic,
                                                       overwrite=overwrite,
                                                     use_multiprocessing=use_multiprocessing,)



        print("len(top_pydcaFNAPC_dict):",len(top_pydcaFNAPC_dict))

        Subject_allPPI_allInfo_dict=dict([((p1,p2),[top_pydcaFNAPC_dict[(p1,p2)][i] for i in topDCA_idx])for p1, p2, len1,len2 in Subject_allPPI_info])
        print("len(Subject_allPPI_allInfo_dict):",len(Subject_allPPI_allInfo_dict))
        Subject_allPPI_allInfo_dict_listDict[(Subject_phylum, Subject_speID)]=Subject_allPPI_allInfo_dict

  

    return(Subject_allPPI_allInfo_dict_listDict)


def get_BestHomologousDCAs_top5DCAs_fromMultiSpes(Query_allPPI_allInfo_dict,
                                    Query_allPPI_top5DCAs_dict,
                                    BestHomologousPP_Subject_allPPI_top5DCAs_listDict,
                                    Query2Subject_BestHomologous_ignoreSubjectDCA_dict_listDict,
                                        topDCA_num=5,
                                    with_status=True,
                                                 ):
    
    Query_BestHomologousDCAs_dict=defaultdict(list)

    for sub_PP in Query_allPPI_allInfo_dict:
        if with_status:
            Query_BestHomologousDCAs_dict[sub_PP]=list(Query_allPPI_allInfo_dict[sub_PP][0])
        else:
            Query_BestHomologousDCAs_dict[sub_PP]=[]
            
        Query_BestHomologousDCAs_dict[sub_PP].extend(Query_allPPI_top5DCAs_dict[sub_PP][0:topDCA_num])

        for Subject_phylum, Subject_speID in Query2Subject_BestHomologous_ignoreSubjectDCA_dict_listDict.keys():
            BestHomologousPP_Subject_allPPI_top5DCAs_dict=BestHomologousPP_Subject_allPPI_top5DCAs_listDict[(Subject_phylum, Subject_speID)]
            Query2Subject_BestHomologous_ignoreSubjectDCA_dict=Query2Subject_BestHomologous_ignoreSubjectDCA_dict_listDict[(Subject_phylum, Subject_speID)]

            best_Subject_PP=Query2Subject_BestHomologous_ignoreSubjectDCA_dict[sub_PP] if sub_PP in Query2Subject_BestHomologous_ignoreSubjectDCA_dict else []
            if len(best_Subject_PP)==0:  # here Subject_PP is tuple with size 2
                Query_BestHomologousDCAs_dict[sub_PP].extend([np.nan for i in range(topDCA_num)])
            else:

                if best_Subject_PP in BestHomologousPP_Subject_allPPI_top5DCAs_dict:
                    Query_BestHomologousDCAs_dict[sub_PP].extend(BestHomologousPP_Subject_allPPI_top5DCAs_dict[best_Subject_PP][0:topDCA_num])

                else:
                    Query_BestHomologousDCAs_dict[sub_PP].extend([np.nan for i in range(topDCA_num)])

                        
    return(Query_BestHomologousDCAs_dict)

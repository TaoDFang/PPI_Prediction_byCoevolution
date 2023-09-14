import argparse
import os 
import glob
import pickle
import multiprocessing as mp
from collections import defaultdict
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd

from IntergrateBestHomologousPPCoEvo_unNameSorted import allQuery2SubjectSingleProteinMapping_getBlasp
from IntergrateBestHomologousPPCoEvo_unNameSorted import chooseBestHomologousPP
from IntergrateBestHomologousPPCoEvo_unNameSorted import Query2Subject_AAmapDic


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()


    parser.add_argument('-q','--Query_tuple', type=str, help='Query_tuple')
    parser.add_argument('-s','--Subject_tupleList', type=str, help='Subject_tupleList')
    parser.add_argument('-seqM','--homologous_SeqMappingPath', type=str, help='homologous_SeqMappingPath')
    parser.add_argument('-m','--homologous_allQuery2SubjectPPIMapping_path', type=str, help='homologous_allQuery2SubjectPPIMapping_path')
    parser.add_argument('-ms','--homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path', type=str, help='homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path')
    parser.add_argument('-mb','--homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path', type=str, help='homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path')


    args = parser.parse_args()
    Query_tuple=args.Query_tuple
    Subject_tupleList=args.Subject_tupleList
    homologous_SeqMappingPath=args.homologous_SeqMappingPath
    homologous_allQuery2SubjectPPIMapping_path=args.homologous_allQuery2SubjectPPIMapping_path
    homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path=args.homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path
    homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path=args.homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path


    Query_tuple=Query_tuple.split("_")
    Query_speID=Query_tuple[1]
    
    Subject_tupleList=Subject_tupleList.split("_")
    Subject_tupleList=[(Subject_tupleList[i],Subject_tupleList[i+1]) for i in range(0,len(Subject_tupleList)-1,2)]
    
    print("Query_tuple:",Query_tuple)
    print("Subject_tupleList:",Subject_tupleList)
    

    with open(homologous_allQuery2SubjectPPIMapping_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_dict_listDict.pickle", 'rb') as handle:
            Query2Subject_QueSpeAllPPI_homologous_dict_listDict=pickle.load(handle)
            
            
    
    with open(homologous_allQuery2SubjectPPIMapping_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict.pickle", 'rb') as handle:
            Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict=pickle.load(handle)

    
    
    # get best  mapping score betwen single proteins 
    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_listDict,Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_blastp_listDict=allQuery2SubjectSingleProteinMapping_getBlasp(Query_speID,Subject_tupleList,homologous_SeqMappingPath,Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict)

    with open(homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_listDict.pickle", 'wb') as handle:
            pickle.dump(Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_listDict,handle)
    with open(homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_blastp_listDict.pickle", 'wb') as handle:
        pickle.dump(Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_blastp_listDict,handle)



#     # AA level mapping.
#     see http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_homologousPPDetection_STRINGPhyBalancePhyla.ipynb
    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Query2Subject_AAmapDic_listDict,Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict=Query2Subject_AAmapDic(Query_speID,Subject_tupleList,
                           homologous_SeqMappingPath,
                           Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict)
    
    with open(homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Query2Subject_AAmapDic_listDict.pickle", 'wb') as handle:
        pickle.dump(Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Query2Subject_AAmapDic_listDict,handle)
    with open(homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict.pickle", 'wb') as handle:
        pickle.dump(Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict,handle)

        
        
    # get best homologous pp mapping 
    Query2Subject_QueSpeAllPPI_BestHomologous_listDict=defaultdict(dict)
    for phylum_speID in Subject_tupleList:
        print(phylum_speID)
        Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_dict=Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_listDict[phylum_speID]


        for Query_PP in Query2Subject_QueSpeAllPPI_homologous_dict_listDict[phylum_speID].keys():

            Subject_PPs=Query2Subject_QueSpeAllPPI_homologous_dict_listDict[phylum_speID][Query_PP]
            Subject_PP=chooseBestHomologousPP(Query_PP,Subject_PPs,Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_dict)

            if len(Subject_PP)==2: # Subject PP is a tuple with two protein so lenght is 2 not 1
                Query2Subject_QueSpeAllPPI_BestHomologous_listDict[phylum_speID][Query_PP]=Subject_PP

    with open(homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_BestHomologous_listDict.pickle", 'wb') as handle:
            pickle.dump(Query2Subject_QueSpeAllPPI_BestHomologous_listDict,handle)    
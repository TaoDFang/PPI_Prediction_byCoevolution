import argparse
import os 
import glob
import pickle
import re
import multiprocessing as mp

from seq_mapping import twoSepMapping


        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('-q','--Query_tuple', type=str, help='Query_tuple')
    parser.add_argument('-s','--Subject_tupleList', type=str, help='Subject_tupleList')
    parser.add_argument('-qb','--QueryProSeqPath_ByProteins', type=str, help='QueryProSeqPath_ByProteins')
    parser.add_argument('-sb','--SubjectProSeqPath_ByProteins', type=str, help='SubjectProSeqPath_ByProteins')
    # parser.add_argument('-mf','--SubjectSpe_MiddleData_folder', type=str, help='SubjectSpe_MiddleData_folder')
    # parser.add_argument('-seqM','--homologous_SeqMappingPath', type=str, help='homologous_SeqMappingPath')
    parser.add_argument('-seqM','--current_homologous_SeqMappingPath', type=str, help='current_homologous_SeqMappingPath')
    parser.add_argument('-m','--homologous_allQuery2SubjectPPIMapping_path', type=str, help='homologous_allQuery2SubjectPPIMapping_path')  
    parser.add_argument('-bp','--blastp_path', type=str, help='blastp_path')  
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')


    args = parser.parse_args()
    Query_tuple=args.Query_tuple
    Subject_tupleList=args.Subject_tupleList
    QueryProSeqPath_ByProteins=args.QueryProSeqPath_ByProteins
    SubjectProSeqPath_ByProteins=args.SubjectProSeqPath_ByProteins
    current_homologous_SeqMappingPath=args.current_homologous_SeqMappingPath
    homologous_allQuery2SubjectPPIMapping_path=args.homologous_allQuery2SubjectPPIMapping_path
    blastp_path=args.blastp_path
    mp_task_nums=int(args.mp_task_nums)
    
    overlap_method="remove"

#     print(SubjectSpe_MiddleData_folder) 
#     # sre=re.compile("([\d+])[^0-9]+")
#     sre=re.compile("([0-9]+)[^0-9]+")
#     Subject_speID,Subject_phylum=sre.findall(SubjectSpe_MiddleData_folder)
#     print("Subject_speID,Subject_phylum:",Subject_speID,Subject_phylum)



    print("QueryProSeqPath_ByProteins:",QueryProSeqPath_ByProteins)
    Query_tuple=Query_tuple.split("_")
    Query_speID=Query_tuple[1]
    

    Subject_tupleList=Subject_tupleList.split("_")
    # Subject_tupleList=[(Subject_tupleList[i],Subject_tupleList[i+1]) for i in range(0,len(Subject_tupleList)-1,2)]
    Subject_spe2phy_dict={Subject_tupleList[i+1]:Subject_tupleList[i] for i in range(0,len(Subject_tupleList)-1,2)}
    
    #before using each keyword, the path  offSubjectProSeqPath_ByProteinSubjectProSeqPath_ByProtein: 411476ByProteins/
    #when use each key workd for , we got full path SubjectProSeqPath_ByProtein: /mnt/mnemo6/tao/nextflow/PPI_Coevolution/STRING_data_11.5/411476ByProteins/
    print("SubjectProSeqPath_ByProteins:",SubjectProSeqPath_ByProteins)
    SubjectProSeqPath_ByProteins_base=os.path.basename(os.path.normpath(SubjectProSeqPath_ByProteins))
    SubjectProSeqPath_ByProteins_base=SubjectProSeqPath_ByProteins_base+"/"
    print("SubjectProSeqPath_ByProteins_base:",SubjectProSeqPath_ByProteins_base)
    sre=re.compile("([0-9]+)[^0-9]+")
    Subject_speID=sre.findall(SubjectProSeqPath_ByProteins_base)[0]
    print("Query_speID,Subject_speID:",Query_speID,Subject_speID)
    Subject_phylum=Subject_spe2phy_dict[Subject_speID]
    print("Subject_phylum:",Subject_phylum)
    
    
    

    
    with open(homologous_allQuery2SubjectPPIMapping_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict.pickle", 'rb') as handle:
        Query2Subject_SubSpeAllPPI_homologous_singleProteinMaping_listDict=pickle.load(handle)
        


    # here there is a problem, new STRIN_rootfodler is a big folder , how to use it in this process witout danger to be trigerd later by other process 
    # unless all new process dont go to the newstrin foot foler anymore 
    # ah ah, use path with wildld card "*", but then the proble is  how to performance for loop ,in nextflow , dont need for loop:
    # check input files secction in https://carpentries-incubator.github.io/workflows-nextflow/05-processes-part1/index.html, use "each" ?

    # old code at : http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_homologousPPDetection.ipynb
    # for phylum, Subject_speID in Subject_tupleList:
    #     print(phylum, Subject_speID)
    #     QueryProSeqPath_ByProteins=newSTRING_rootFolder+Query_speID+"ByProteins/"
    #     SubjectProSeqPath_ByProteins=newSTRING_rootFolder+Subject_speID+"ByProteins/"

#     current_homologous_SeqMappingPath=homologous_SeqMappingPath+"EggNogMaxLevel2_QuerySpe_ID"+Query_speID+"and"+"SubjectSpe_ID"+Subject_speID+'/'
#     print(current_homologous_SeqMappingPath)

#     if not os.path.exists(current_homologous_SeqMappingPath):
#         os.makedirs(current_homologous_SeqMappingPath)

    current_homologous_SeqMappingPath_results=glob.glob(current_homologous_SeqMappingPath+"*.pickle")
    current_homologous_SeqMappingPath_resultsDic={f:1 for f in current_homologous_SeqMappingPath_results}


    temp_Query2Subject_SubSpeAllPPI_homologous_singleProteinMaping_dict=Query2Subject_SubSpeAllPPI_homologous_singleProteinMaping_listDict[(Subject_phylum,Subject_speID)]

    ArgFor_twoSepMapping=[(Query_proID,Subject_proID,overlap_method,QueryProSeqPath_ByProteins, SubjectProSeqPath_ByProteins,current_homologous_SeqMappingPath,blastp_path) for Query_proID,Subject_proIDs in temp_Query2Subject_SubSpeAllPPI_homologous_singleProteinMaping_dict.items() for Subject_proID in Subject_proIDs]
    print("len(ArgFor_twoSepMapping):",len(ArgFor_twoSepMapping))

    ArgFor_twoSepMapping=[(Query_proID,Subject_proID,overlap_method,Query_fastaPath, Subject_fastaPath,outputPath,blastp_path) for Query_proID,Subject_proID,overlap_method,Query_fastaPath, Subject_fastaPath,outputPath,blastp_path in ArgFor_twoSepMapping if (outputPath+Query_proID+"and"+Subject_proID+"Subject2Query_AAmapDic_overlapMethod"+overlap_method+".pickle") not in current_homologous_SeqMappingPath_resultsDic ]
    print("len(ArgFor_twoSepMapping):",len(ArgFor_twoSepMapping))

    pool=mp.Pool(mp_task_nums)
    pool.map(twoSepMapping,ArgFor_twoSepMapping)
    pool.close()
    # time.sleep(1)


        
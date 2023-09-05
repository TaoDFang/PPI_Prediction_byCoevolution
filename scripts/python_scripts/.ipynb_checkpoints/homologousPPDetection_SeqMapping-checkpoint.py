import argparse
import os 
import glob
import pickle
import multiprocessing as mp

from seq_mapping import twoSepMapping

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()


    parser.add_argument('-seqM','--homologous_SeqMappingPath', type=str, help='homologous_SeqMappingPath')


    args = parser.parse_args()
    homologous_SeqMappingPath=args.homologous_SeqMappingPath





#     Query_speID=Query_tuple[1]
#     overlap_method="remove"



#     # here there is a problem, new STRIN_rootfodler is a big folder , how to use it in this process witout danger to be trigerd later by other process 
#     # unless all new process dont go to the newstrin foot foler anymore 
#     # ah ah, use path with wildld card "*", but then the proble is  how to performance for loop ,in nextflow , dont need for loop:
#     # check input files secction in https://carpentries-incubator.github.io/workflows-nextflow/05-processes-part1/index.html ,use "each" ?

#     # old code at : http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_homologousPPDetection.ipynb
#     # for phylum, Subject_speID in Subject_tupleList:
#     #     print(phylum, Subject_speID)
#     #     QueryProSeqPath_ByProteins=newSTRING_rootFolder+Query_speID+"ByProteins/"
#     #     SubjectProSeqPath_ByProteins=newSTRING_rootFolder+Subject_speID+"ByProteins/"

#         current_homologous_SeqMappingPath=homologous_SeqMappingPath+"EggNogMaxLevel2_QuerySpe_ID"+Query_speID+"and"+"SubjectSpe_ID"+Subject_speID+'/'
#         print(current_homologous_SeqMappingPath)

#         current_homologous_SeqMappingPath_results=glob.glob(current_homologous_SeqMappingPath+"*.pickle")
#         current_homologous_SeqMappingPath_resultsDic={f:1 for f in current_homologous_SeqMappingPath_results}


#         if not os.path.exists(current_homologous_SeqMappingPath):
#             os.makedirs(current_homologous_SeqMappingPath)

#         temp_Query2Subject_SubSpeAllPPI_homologous_singleProteinMaping_dict=Query2Subject_SubSpeAllPPI_homologous_singleProteinMaping_listDict[(phylum,Subject_speID)]

#         ArgFor_twoSepMapping=[(Query_proID,Subject_proID,overlap_method,QueryProSeqPath_ByProteins, SubjectProSeqPath_ByProteins,current_homologous_SeqMappingPath) for Query_proID,Subject_proIDs in temp_Query2Subject_SubSpeAllPPI_homologous_singleProteinMaping_dict.items() for Subject_proID in Subject_proIDs]
#         print("len(ArgFor_twoSepMapping):",len(ArgFor_twoSepMapping))

#         ArgFor_twoSepMapping=[(Query_proID,Subject_proID,overlap_method,Query_fastaPath, Subject_fastaPath,outputPath) for Query_proID,Subject_proID,overlap_method,Query_fastaPath, Subject_fastaPath,outputPath in ArgFor_twoSepMapping if (outputPath+Query_proID+"and"+Subject_proID+"Subject2Query_AAmapDic_overlapMethod"+overlap_method+".pickle") not in current_homologous_SeqMappingPath_resultsDic ]
#         print("len(ArgFor_twoSepMapping):",len(ArgFor_twoSepMapping))

#         pool=mp.Pool(20)
#         pool.map(twoSepMapping,ArgFor_twoSepMapping)
#         pool.close()
#         time.sleep(1)


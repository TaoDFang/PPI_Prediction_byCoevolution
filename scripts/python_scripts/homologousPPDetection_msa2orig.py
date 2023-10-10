import argparse
import os 

import pickle



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()


    parser.add_argument('-q','--Query_tuple', type=str, help='Query_tuple')
    parser.add_argument('-s','--Subject_tupleList', type=str, help='Subject_tupleList')
    parser.add_argument('-ms','--homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path', type=str, help='homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path')


    args = parser.parse_args()
    Query_tuple=args.Query_tuple
    Subject_tupleList=args.Subject_tupleList


        
    with open(homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict.pickle", 'rb') as handle:
        Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict=pickle.dump(handle)

    Query_allProteins=set()
    Subject_allProteins=defaultdict(set)
    for phylum_speID,temp_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict in Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict.items():
        for Query_proID,Subject_proID in temp_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict:
            Query_allProteins.add(Query_proID)
            Subject_allProteins[phylum_speID].add(Subject_proID)

    print(len(Query_allProteins))

    for phylum_speID,currentSubject_proteins in Subject_allProteins.items():
        print(phylum_speID)
        print(len(currentSubject_proteins))


    Query_allProteins_msa2orig_dict=dict()
    for Query_pro in Query_allProteins:
        Query_pro_trackGapsPos_file=Query_msa_trackGapsPos_path+Query_pro+'_allSeq_keptAA_posInOrigalPro.pickle'
        if os.path.exists(Query_pro_trackGapsPos_file): # not all protein passed filtering step to final paired MSA 
            with open(Query_pro_trackGapsPos_file, 'rb') as handle:
                protein_posTrack=pickle.load(handle)
            protein_origPos=protein_posTrack[0][1]
            Query_allProteins_msa2orig_dict[Query_pro]=dict(zip(range(len(protein_origPos)),protein_origPos))

        else:
            Query_allProteins_msa2orig_dict[Query_pro]=None


    Subject_allProteins_msa2orig_listDict=defaultdict(dict)
    for phylum_speID,currentSubject_proteins in Subject_allProteins.items():
        print(phylum_speID)
        currentSubject_msa_trackGapsPos_path=newSTRING_rootFolder+phylum_speID[1]+"_EggNOGmaxLevel"+phylum_speID[0]+"_newSingleMSA_hmmalign_removeGaps_trackGapPos/"   
        for Subject_pro in currentSubject_proteins:
            Subject_pro_trackGapsPos_file=currentSubject_msa_trackGapsPos_path+Subject_pro+'_allSeq_keptAA_posInOrigalPro.pickle'
            if os.path.exists(Subject_pro_trackGapsPos_file):
                with open(Subject_pro_trackGapsPos_file, 'rb') as handle:
                    protein_posTrack=pickle.load(handle)
                protein_origPos=protein_posTrack[0][1]
                Subject_allProteins_msa2orig_listDict[phylum_speID][Subject_pro]=dict(zip(range(len(protein_origPos)),protein_origPos))
            else:
                Subject_allProteins_msa2orig_listDict[phylum_speID][Subject_pro]=None


    with open(homologous_msa2origMappingPath+"QueSpeAllPPI"+Query_tuple[1]+"_allProteins_msa2orig_dict.pickle", 'wb') as handle:
            pickle.dump(Query_allProteins_msa2orig_dict,handle)


    with open(homologous_msa2origMappingPath+"QueSpeAllPPI"+"Subject_allProteins_msa2orig_listDict.pickle", 'wb') as handle:
            pickle.dump(Subject_allProteins_msa2orig_listDict,handle)
        



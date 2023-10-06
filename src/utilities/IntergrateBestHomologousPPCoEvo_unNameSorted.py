import os 
import pickle


from collections import defaultdict
from Bio.Blast import NCBIXML
import numpy as np
import pandas as pd

def getEggNOG_homologousNameUnSortedPP_COG2PP_reversedPP_dict(record):
#                                                  QuerySpe_ID="511145",
#                                                  QuerySpe_ID=None,
#                                                  max_level="2",
#                                                  EggNOG_groupPath="/mnt/mnemo6/tao/STRING_Data_11.5/eggnog5AddSTRING11.5_Species/groups/",
#                                                  homologous_COG2PP_path="/mnt/mnemo6/tao/PPI_Coevolution/STRING_data_11.5/homologous_pp/"

    CurrentSpe_ID,max_level,EggNOG_groupPath,homologous_COG2PP_path=record
    
    EggNOG_homologousPP_COG2PP_CurrentSpe_dict_file=homologous_COG2PP_path+"EggNogMaxLevel"+max_level+"_"+CurrentSpe_ID+'_HomologousCOG2PP_NameUnSorted_dict.pickle'
    
    if not os.path.exists(EggNOG_homologousPP_COG2PP_CurrentSpe_dict_file): # this is not working in nextflow but useful when run seperately 
        EggNOG_group_level2=pd.read_csv(EggNOG_groupPath+max_level+".tsv",
                                        header=None,index_col=None,sep="\t")

        EggNOG_group_level2_speID=EggNOG_group_level2.iloc[:,2].values.tolist()
        EggNOG_group_level2_speID=[str(sid) for sid in EggNOG_group_level2_speID]

        EggNOG_group_level2_idx=[idx for idx,sid in  enumerate(EggNOG_group_level2_speID) if (sid==CurrentSpe_ID)]
        EggNOG_group_level2=EggNOG_group_level2.iloc[EggNOG_group_level2_idx,:]
        EggNOG_group_level2=EggNOG_group_level2.sort_values(by=1)




        # create a dictiontionay whose key is eggnog group id/name and keys are speceis and protein belong to this eggnog group 
        EggNOG_group_level2_groupDict=defaultdict(list)

        EggNOG_group_level2_list=EggNOG_group_level2.values.tolist()

        for _,group_id,sid,pname in EggNOG_group_level2_list:
            EggNOG_group_level2_groupDict[group_id].append((str(sid),pname))

        print(CurrentSpe_ID,"EggNOG_group_level2_groupDict",len(EggNOG_group_level2_groupDict))

        # for all  eggnog group pairs , get homologous protein pairs belong to these eggnog group pairs 
        # and creat two dictinoary to store such information for query and Query species 
        
        # here have to order eggnog group. so each eggnog group pair, only one recored remained 
        EggNOG_group_level2_uniqueGroups=sorted(list(EggNOG_group_level2_groupDict.keys()))
        EggNOG_homologousPP_COG2PP_CurrentSpe_dict=defaultdict(list)

        for ig in range(len(EggNOG_group_level2_uniqueGroups)-1):
            ig_eggnogID=EggNOG_group_level2_uniqueGroups[ig]
            ig_values = EggNOG_group_level2_groupDict[ig_eggnogID]
            for jg  in range(ig+1,len(EggNOG_group_level2_uniqueGroups)):
                jg_eggnogID=EggNOG_group_level2_uniqueGroups[jg]
                jg_values = EggNOG_group_level2_groupDict[jg_eggnogID]
                for ig_sid, ig_pname in ig_values:
                    for jg_sid,jg_pname in jg_values:
                        #if (ig_sid == jg_sid):  # actually this step is not neceeeay as of course they are from all species
                        # because we order COG groups before 
                        # there we make sure both (protein1, protein2 ) and (protein2, protein1 )  have records
                        EggNOG_homologousPP_COG2PP_CurrentSpe_dict[(ig_eggnogID,jg_eggnogID)].append((ig_pname,jg_pname))
                        EggNOG_homologousPP_COG2PP_CurrentSpe_dict[(jg_eggnogID,ig_eggnogID)].append((jg_pname,ig_pname))
        print("len(EggNOG_homologousPP_COG2PP_CurrentSpe_dict:",len(EggNOG_homologousPP_COG2PP_CurrentSpe_dict))

        

        #notice even for nextlfow , the folder homologous_COG2PP_path is temperory, it need to be creaed before to avoid file not exist problem
        with open(EggNOG_homologousPP_COG2PP_CurrentSpe_dict_file, 'wb') as handle:
            pickle.dump(EggNOG_homologousPP_COG2PP_CurrentSpe_dict, handle)


            

def allQuery2SubjectPPIMapping_getNameUnSorted_COGasLink(Query_tuple,Subject_tupleList,PPIInfoBeforeCoEvoComp_csv,homologous_COG2PP_path):
    
    # use  postive and negative ppi to only use pp that have large Nf90 value, and after removing deep homologs  
    currentSpe_allPPIs_beforeCoEvoComp_frame = pd.read_csv(PPIInfoBeforeCoEvoComp_csv,header=0,index_col=None,sep="\t")
    print("currentSpe_allPPIs_beforeCoEvoComp_frame.shape:",currentSpe_allPPIs_beforeCoEvoComp_frame.shape)
    currentSpe_allPPIs_beforeCoEvoComp_info=currentSpe_allPPIs_beforeCoEvoComp_frame.values.tolist()
    allPPIs=[(p1,p2) for p1, p2,_ in currentSpe_allPPIs_beforeCoEvoComp_info]
    print("len(allPPIs):",len(allPPIs))
    Query_allPPI_allInfo_dict={pp:1 for pp in allPPIs}
    print("len(Query_allPPI_allInfo_dict)",len(Query_allPPI_allInfo_dict))
    
    
   # this is cog to pp relationship in eggnot level 2 despite spceis
    Query_EggNOG_homologousPP_COG2PP_dict_file=homologous_COG2PP_path+"EggNogMaxLevel2_"+Query_tuple[1]+'_HomologousCOG2PP_NameUnSorted_dict.pickle'
    with open(Query_EggNOG_homologousPP_COG2PP_dict_file, 'rb') as handle:
            Query_EggNOG_homologousPP_COG2PP_dict=pickle.load(handle)    
    print("len(Query_EggNOG_homologousPP_COG2PP_dict):",len(Query_EggNOG_homologousPP_COG2PP_dict))


    Query_EggNOG_homologousPP_PP2COG_dict=defaultdict(list)
    for k, vs in Query_EggNOG_homologousPP_COG2PP_dict.items():
        for v in vs:
            Query_EggNOG_homologousPP_PP2COG_dict[v].append(k)
    print(len(Query_EggNOG_homologousPP_PP2COG_dict))

    #keep in mind here "PP2COG_dict" is already contain reversed records
    # here filter cog to pp relation by only keep pp with calculated high dca scores to save computational time later 
    Query_EggNOG_homologousPP_PP2COG_dict={ pp: Query_EggNOG_homologousPP_PP2COG_dict[pp] for pp in Query_allPPI_allInfo_dict}
    print("len(Query_EggNOG_homologousPP_PP2COG_dict):",len(Query_EggNOG_homologousPP_PP2COG_dict))

    
    # for subject data 
    Query2Subject_homologous_dict_listDict=dict()
    for Subject_phylum, Subject_speID in Subject_tupleList:
        print("Subject_phylum, Subject_speID:",Subject_phylum, Subject_speID)

        # this is cog to pp relationship in eggnot level 2 despite spceis
        Subject_EggNOG_homologousPP_COG2PP_dict_file=homologous_COG2PP_path+"EggNogMaxLevel2_"+Subject_speID+'_HomologousCOG2PP_NameUnSorted_dict.pickle'
        with open(Subject_EggNOG_homologousPP_COG2PP_dict_file, 'rb') as handle:
                Subject_EggNOG_homologousPP_COG2PP_dict=pickle.load(handle)
        print(len(Subject_EggNOG_homologousPP_COG2PP_dict))



        # combined info
        Query2Subject_homologous_dict=defaultdict(list)
        for Query_pp, cogs in Query_EggNOG_homologousPP_PP2COG_dict.items():
            for cog in cogs: 
                if cog in Subject_EggNOG_homologousPP_COG2PP_dict:
                    Subject_pps=Subject_EggNOG_homologousPP_COG2PP_dict[cog]
                    # for each cog , should be a append 
                    Query2Subject_homologous_dict[Query_pp].extend(Subject_pps)
                    #Namesorted_Query2Subject_homologous_dict[subject_pp]=query_pps

        #thats because in query speceis, some different  protein pairs can exist  in diffrent cog pairs 
        #that can be mapped to same protein pairs in subject species 
        # so bacically , one protein pair in sub species and one protein pair in query sepceis share several cog paris 
        print(len(Query2Subject_homologous_dict))
        Query2Subject_homologous_dict={k: list(set(v)) for k, v in Query2Subject_homologous_dict.items()}
        
        print(len(Query2Subject_homologous_dict))

        Query2Subject_homologous_dict_listDict[(Subject_phylum, Subject_speID)]=Query2Subject_homologous_dict
    


    return(Query2Subject_homologous_dict_listDict)


def allQuery2SubjectSingleProteinMapping(Query2Subject_QueSpeAllPPI_homologous_dict_listDict):
 

    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict=defaultdict(dict)
    for phylum_speID,temp_Query2Subject_QueSpeAllPPI_homologous_dict in Query2Subject_QueSpeAllPPI_homologous_dict_listDict.items():
        print(phylum_speID)
        for Query_pp, Subject_pps in temp_Query2Subject_QueSpeAllPPI_homologous_dict.items():
            Query_pp1,Query_pp2=Query_pp
            for Subject_pp in Subject_pps:
    #             if(len(Subject_pps)==1):
    #                 print(Subject_pp)
                Subject_pp1,Subject_pp2=Subject_pp

                # here we should notice one  single protein can have multiple homologous protein pairs
                # and same one singel protein to single protein mapping can appear many times in different homologous protein pair mapping 
                if Query_pp1 not in Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID]:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID][Query_pp1]={Subject_pp1} #set
                else:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID][Query_pp1].add(Subject_pp1)

                if Query_pp2 not in Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID]:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID][Query_pp2]={Subject_pp2} # set 
                else:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID][Query_pp2].add(Subject_pp2)


    return(Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict)

def allQuery2SubjectSingleProteinMapping_getBlasp(Query_speID,Subject_tupleList,homologous_SeqMappingPath,Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict):    
    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_listDict=defaultdict(dict)
    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_blastp_listDict=defaultdict(dict)
    for phylum, Subject_speID in Subject_tupleList:
        print(phylum, Subject_speID)

        current_homologous_SeqMappingPath=homologous_SeqMappingPath+"EggNogMaxLevel2_QuerySpe_ID"+Query_speID+"and"+"SubjectSpe_ID"+Subject_speID+'/'

        temp_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_dict=Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[(phylum,Subject_speID)]

        for Query_proID,Subject_proIDs in temp_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_dict.items():
            for Subject_proID in Subject_proIDs:
                #print(current_homologous_SeqMappingPath+Query_proID+"and"+Subject_proID+".xml")
                with open(current_homologous_SeqMappingPath+Query_proID+"and"+Subject_proID+".xml") as handle:
                    blast_record = NCBIXML.read(handle)
                    #print( blast_record.num_hits)

                    # ???? here there is a problem is that if blast_record.alignments is None
                    # then this for loop will never get started and thus sorted_alignment_hsps and beste_hsp
                    # will just be the values from last for loop 
                    #print(len(blast_record.alignments))
                    if len(blast_record.alignments)>0:
                        for alignment in blast_record.alignments:
                            #print("times")
                            sorted_alignment_hsps=sorted(alignment.hsps,key=lambda hsp: hsp.bits,reverse=True) 
                            best_hsp=sorted_alignment_hsps[0]
                            #print("best_hsp.bits:",best_hsp.bits)
                        Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_listDict[(phylum,Subject_speID)][(Query_proID,Subject_proID)]=best_hsp.bits
                    else:
                        #print("here empty blastp results")
                        Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_listDict[(phylum,Subject_speID)][(Query_proID,Subject_proID)]=np.nan

                # '6 qseqid  qlen sseqid  slen qstart qend sstart send evalue bitscore pident qcovs qcovhsp'
                try:
                    out_frame=pd.read_csv(current_homologous_SeqMappingPath+Query_proID+"and"+Subject_proID+".txt",header=None,sep="\t")
                    out_frame=out_frame.sort_values(by=[9],ascending=True) # sorted by bit scores 
                    out_values=out_frame.values.tolist()
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_blastp_listDict[(phylum,Subject_speID)][(Query_proID,Subject_proID)]=out_values

                except:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_blastp_listDict[(phylum,Subject_speID)][(Query_proID,Subject_proID)]=np.nan
                    
    return (Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_bestHspBits_listDict,
           Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_blastp_listDict)



def Query2Subject_AAmapDic(Query_speID,Subject_tupleList,
                           homologous_SeqMappingPath,
                           Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict,
                           overlap_method="remove"):


    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Query2Subject_AAmapDic_listDict=defaultdict(dict)
    for phylum, Subject_speID in Subject_tupleList:

        current_homologous_SeqMappingPath=homologous_SeqMappingPath+"EggNogMaxLevel2_QuerySpe_ID"+Query_speID+"and"+"SubjectSpe_ID"+Subject_speID+'/'


        temp_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_dict=Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[(phylum,Subject_speID)]

        for Query_proID,Subject_proIDs in temp_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_dict.items():
            for Subject_proID in Subject_proIDs:
                with open(current_homologous_SeqMappingPath+Query_proID+"and"+Subject_proID+"subject2query_AAmapDic_overlapMethod"+overlap_method+".pickle", 'rb') as handle:
                    Subject2Query_AAmapDic=pickle.load(handle)

                    # here pos in AAmapDic start with 1
                    # !! we force it to start from 0 , so its consistent with  0-starting pos from other object
                    Subject2Query_AAmapDic={k-1:v-1for k, v in Subject2Query_AAmapDic.items()}

                    Query2Subject_AAmapDic={v:k for k , v in Subject2Query_AAmapDic.items()}

                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Query2Subject_AAmapDic_listDict[(phylum,Subject_speID)][(Query_proID,Subject_proID)]=Query2Subject_AAmapDic


    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict=defaultdict(dict)
    for phylum, Subject_speID in Subject_tupleList:

        current_homologous_SeqMappingPath=homologous_SeqMappingPath+"EggNogMaxLevel2_QuerySpe_ID"+Query_speID+"and"+"SubjectSpe_ID"+Subject_speID+'/'


        temp_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_dict=Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[(phylum,Subject_speID)]

        for Query_proID,Subject_proIDs in temp_Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_dict.items():
            for Subject_proID in Subject_proIDs:
                with open(current_homologous_SeqMappingPath+Query_proID+"and"+Subject_proID+"subject2query_AAmapDic_overlapMethod"+overlap_method+".pickle", 'rb') as handle:
                    Subject2Query_AAmapDic=pickle.load(handle)

                    # here pos in AAmapDic start with 1
                    # !! we force it to start from 0 , so its consistent with  0-starting pos from other object
                    Subject2Query_AAmapDic={k-1:v-1for k, v in Subject2Query_AAmapDic.items()}

                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict[(phylum,Subject_speID)][(Query_proID,Subject_proID)]=Subject2Query_AAmapDic


    return(Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Query2Subject_AAmapDic_listDict,Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_Subject2Query_AAmapDic_listDict)


def chooseBestHomologousPP(Query_PP,Subject_PPs,Query2Subject_homologous_ignoreSubjectDCA_singleProteinMaping_Subject2Query_bestHspBits_dict):
    
    #print(Subject_PPs)
    BestHomologousPP=()
    BestHomologousPP_bestHspBits=()
    # when Subject_PPs is none, for lood dont stat and empty Subject pp returned 
    for Subject_PP in Subject_PPs:
        cur_QuerySubjectPP1_bestHspBits=Query2Subject_homologous_ignoreSubjectDCA_singleProteinMaping_Subject2Query_bestHspBits_dict[(Query_PP[0],Subject_PP[0])]
        cur_QuerySubjectPP2_bestHspBits=Query2Subject_homologous_ignoreSubjectDCA_singleProteinMaping_Subject2Query_bestHspBits_dict[(Query_PP[1],Subject_PP[1])]

        if not ((np.isnan(cur_QuerySubjectPP1_bestHspBits)) or (np.isnan(cur_QuerySubjectPP2_bestHspBits))) :
            #print(Subject_PP)
            #print(cur_QuerySubjectPP1_bestHspBits)
            #print(cur_QuerySubjectPP2_bestHspBits)
            
            if len(BestHomologousPP)==0:
                BestHomologousPP=Subject_PP
                BestHomologousPP_bestHspBits=(cur_QuerySubjectPP1_bestHspBits,cur_QuerySubjectPP2_bestHspBits)
            else:
                # compare lower bestHspBits of currrent Subject pp and homologous pp
                # choose one with high lower bestHspbit
                if min(BestHomologousPP_bestHspBits)< min(cur_QuerySubjectPP1_bestHspBits,cur_QuerySubjectPP2_bestHspBits):
                    BestHomologousPP=Subject_PP
                    BestHomologousPP_bestHspBits=(cur_QuerySubjectPP1_bestHspBits,cur_QuerySubjectPP2_bestHspBits)
                # when lower bestHspBits  are same , compare another protein 
                elif (min(BestHomologousPP_bestHspBits)== min(cur_QuerySubjectPP1_bestHspBits,cur_QuerySubjectPP2_bestHspBits)) and (max(BestHomologousPP_bestHspBits)< max(cur_QuerySubjectPP1_bestHspBits,cur_QuerySubjectPP2_bestHspBits)):
                    BestHomologousPP=Subject_PP
                    BestHomologousPP_bestHspBits=(cur_QuerySubjectPP1_bestHspBits,cur_QuerySubjectPP2_bestHspBits) 
                    
    return(BestHomologousPP)





def getMetaFrame_FullBestHomologousPP(Query_tuple,EggNOG_maxLevel,currentSpe_TaxID,STRING_Version,
                                      benchmark_suffix="STRINPhyPPI_Benchmark/",
                                     smallPhylum="",
                                     BestHomologousPP_filePrefix="BestHomologousPP",
                                        CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING11.5/",
                                     ):
    '''
    here full BestHomologousPP means all BestHomologous pp compared with Query speceis
    It does not include all pp in origin indepenpendt benchmark dataset 
    '''
    #CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING"+STRING_Version+"/"
    input_root_folder=CoEvo_data_folder+smallPhylum+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_eggNOGfilteredData/"

    
    Subect_BestHomologousPP_prefix="BestHomologousPPFor"+Query_tuple[1]+"AtEggNOGmaxLevel"+Query_tuple[0]+"_"
    BestHomologousPP_Benchmark_folder=input_root_folder+Subect_BestHomologousPP_prefix+benchmark_suffix #
    print("BestHomologousPP_Benchmark_folder:",BestHomologousPP_Benchmark_folder)
        
    allPPI_allInfo_frame=pd.read_csv(BestHomologousPP_Benchmark_folder+BestHomologousPP_filePrefix+"_allPPI_allInfo_frame.csv",
                                 header=0,index_col=None,sep="\t")
    
    
    print("allPPI_allInfo_frame.shape:",allPPI_allInfo_frame.shape)
    
    if benchmark_suffix == "STRINPhyPPI_Benchmark/":
        Pos_allPPI_allInfo_frame=allPPI_allInfo_frame.loc[allPPI_allInfo_frame['benchmark_status']=="P",:]
        #print("Pos_allPPI_allInfo_frame.shape:",Pos_allPPI_allInfo_frame.shape)
        Pos_allPPI_allInfo_frame=Pos_allPPI_allInfo_frame.sort_values(by="maxBetDCA_score",ascending=False)

        Neg_allPPI_allInfo_frame=allPPI_allInfo_frame.loc[allPPI_allInfo_frame['benchmark_status']=="N",:]
        #print("Neg_allPPI_allInfo_frame.shape:",Neg_allPPI_allInfo_frame.shape)

        Neg_allPPI_allInfo_frame=Neg_allPPI_allInfo_frame.sort_values(by="maxBetDCA_score",ascending=False)

        allPPI_allInfo_frame= pd.concat([Pos_allPPI_allInfo_frame,Neg_allPPI_allInfo_frame])
    elif benchmark_suffix == "KEGG_Benchmark/":
        allPPI_allInfo_frame=allPPI_allInfo_frame.sort_values(by="maxBetDCA_score",ascending=False)
    elif benchmark_suffix == "QianCong_Benchmark/":
        allPPI_allInfo_frame=allPPI_allInfo_frame.sort_values(by="maxBetDCA_score",ascending=False)
        
    else:
        return(allPPI_allInfo_frame)
    

    
    return(allPPI_allInfo_frame)






def collect_BestHomologousDCAs_OneSpeOneScore_OnlyTopPosNeg(allPPI_allInfoFrame,
                                                      Query_BestHomologousDCAs_dict,
                                                      ML_inputPath,CoEvo_type="DCA",
                                                           BestHomologousDCAs_dict_withStatus=True):


    X_fileName=ML_inputPath+"OnlyTopPosNeg_NonPara_X_BestHomologous"+CoEvo_type+"s_OneSpeOneScore.npy"
    Y_fileName=ML_inputPath+"OnlyTopPosNeg_NonPara_Y_BestHomologous"+CoEvo_type+"s_OneSpeOneScore.npy"
    #if not os.path.exists(X_fileName):



    NonPara_allPPI_info=allPPI_allInfoFrame.loc[:,["STRING_ID1","STRING_ID2","benchmark_status"]].values.tolist()



    if BestHomologousDCAs_dict_withStatus:
        OnlyTopPosNeg_NonPara_X=np.zeros((len(NonPara_allPPI_info),len(list(Query_BestHomologousDCAs_dict.values())[0])-1))
        OnlyTopPosNeg_NonPara_Y=np.zeros((len(NonPara_allPPI_info),1))
        for idx,l in enumerate(NonPara_allPPI_info):
            p1, p2 , label=l
            numeric_lable=1 if label=="P" else 0

            OnlyTopPosNeg_NonPara_X[idx,:]=Query_BestHomologousDCAs_dict[(p1,p2)][1:]
            OnlyTopPosNeg_NonPara_Y[idx]=numeric_lable
    else:
        OnlyTopPosNeg_NonPara_X=np.zeros((len(NonPara_allPPI_info),len(list(Query_BestHomologousDCAs_dict.values())[0])))
        OnlyTopPosNeg_NonPara_Y=np.zeros((len(NonPara_allPPI_info),1))
        for idx,l in enumerate(NonPara_allPPI_info):
            p1, p2 , label=l
            numeric_lable=1 if label=="P" else 0

            OnlyTopPosNeg_NonPara_X[idx,:]=Query_BestHomologousDCAs_dict[(p1,p2)]
            OnlyTopPosNeg_NonPara_Y[idx]=numeric_lable


    np.save(X_fileName,OnlyTopPosNeg_NonPara_X)
    np.save(Y_fileName,OnlyTopPosNeg_NonPara_Y)

    OnlyTopPosNeg_NonPara_XtopCoEvos=np.load(X_fileName)

    OnlyTopPosNeg_NonPara_YtopCoEvos=np.load(Y_fileName)
    OnlyTopPosNeg_NonPara_YtopCoEvos=np.reshape(OnlyTopPosNeg_NonPara_YtopCoEvos,(OnlyTopPosNeg_NonPara_YtopCoEvos.shape[0],))#
    
    return(OnlyTopPosNeg_NonPara_XtopCoEvos,OnlyTopPosNeg_NonPara_YtopCoEvos)
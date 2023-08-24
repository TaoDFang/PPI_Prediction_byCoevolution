#!/usr/bin/env python
# coding: utf-8


import argparse
import os 
import pandas as pd

parser = argparse.ArgumentParser(description='RemoveRedundantProteins')
parser.add_argument('-r','--redundant_proteins_csvFile', type=str, help='redundant_proteins_csvFile')
parser.add_argument('-b','--currentSpe_withinBlastPath', type=str, help='currentSpe_withinBlastPath')


args = parser.parse_args()
redundant_proteins_csvFile=args.redundant_proteins_csvFile
currentSpe_withinBlastPath=args.currentSpe_withinBlastPath


if not os.path.exists(redundant_proteins_csvFile):
    withinBlast_result=pd.read_csv(currentSpe_withinBlastPath+"all2all.txt",comment='#',header=None,sep="\t")

    #remove blast reuslt of itself
    withinBlast_result=withinBlast_result.loc[withinBlast_result.iloc[:,0]!=withinBlast_result.iloc[:,3],:]

    withinBlast_result_uniPros=list(set(withinBlast_result.iloc[:,0].values.tolist()))


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
    unique_redudant_tupleList=set(redudant_tupleList)
    print("len(unique_redudant_tupleList):",len(unique_redudant_tupleList))

    redudant_short_pros=[temp[1] for temp in unique_redudant_tupleList] # notice here temp[1] are already shorter redunat protein
    print("len(redudant_short_pros)",len(redudant_short_pros))

    unique_redudant_short_pros=set(redudant_short_pros)
    print("len(unique_redudant_short_pros)",len(unique_redudant_short_pros))


    unique_redudant_short_pros_frame=pd.DataFrame(unique_redudant_short_pros,columns=["pro_name"])
    unique_redudant_short_pros_frame.to_csv(redundant_proteins_csvFile,
                                            header=True,index=None)
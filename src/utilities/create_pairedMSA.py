import pandas as pd 
import numpy as np 
import math
import glob
import pickle

import os 
import sys
import re

from Bio import AlignIO
from Bio import SeqIO


def getComStrSpeForPairedPro_Save(spe_TaxID,hmmalign_path,QueryPro_pro1,QueryPro_pro2,pairedMSA_unfiltered_folder):  
    
    try:
        alignment1 = AlignIO.read(hmmalign_path+QueryPro_pro1+".fa", "fasta")
        alignment2 = AlignIO.read(hmmalign_path+QueryPro_pro2+".fa", "fasta")

        retemp=re.compile("([0-9]*)\..*")

        paired_msa_list=[]


        proNames1=[record.id for record in alignment1]
        StrSpeID1=[int(retemp.match(pro).group(1)) for pro in proNames1]


        proNames2=[record.id for record in alignment2]
        StrSpeID2=[int(retemp.match(pro).group(1)) for pro in proNames2]

        comStrSpe=list(set(StrSpeID1).intersection(set(StrSpeID2)))
        comStrSpe.insert(0,comStrSpe.pop(comStrSpe.index(int(spe_TaxID))))


        if len(comStrSpe)>1:

            align1_dict=dict()
            for record in alignment1:
                align1_dict[record.id]=str(record.seq)

            align2_dict=dict()
            for record in alignment2:
                align2_dict[record.id]=str(record.seq)


            for spe in comStrSpe:
                OG_Pro1=proNames1[StrSpeID1.index(spe)]
                OG_Pro2=proNames2[StrSpeID2.index(spe)]

                paired_msa_list.append([OG_Pro1+"and"+OG_Pro2,align1_dict[OG_Pro1]+align2_dict[OG_Pro2]])




        if len(paired_msa_list)>1:
            OG_output_str=str()
            for name,seq in paired_msa_list:
                OG_output_str=OG_output_str+">"+name+"\n"+str(seq)+"\n"
            with open(pairedMSA_unfiltered_folder+QueryPro_pro1+"and"+QueryPro_pro2+".fasta", "w") as myfile:
                myfile.write(OG_output_str)


            QueryPro1_msa_len,QueryPro2_msa_len=len(align1_dict[OG_Pro1]),len(align2_dict[OG_Pro2])
            fake_Nf90=len(paired_msa_list)/(math.sqrt(QueryPro1_msa_len+QueryPro2_msa_len))

            return((QueryPro1_msa_len,QueryPro2_msa_len,fake_Nf90,paired_msa_list))


        else:
            print("not enough comon species")
            return (None)


    # except OSError as e:
    except Exception  as e:  # By this way we can know about the type of error occurring
        print("weird thing happens here, this two MSA dont contain  current species as thier common species ",QueryPro_pro1,QueryPro_pro2,e)
        # print("len(alignment1),len(alignment2):",len(alignment1),len(alignment2))
        #happend when single msa of  one of protein is missing (didnt pass filtering steps)
    

        
        
def HHfilter2Save2Transfer2Fasta2Save(QueryPro_pro1,QueryPro_pro2,
                                      QueryPro1_msa_len,QueryPro2_msa_len,
                                      paired_msa,
                                      Nf90_Thres,
                                      pairedMSA_unfiltered_folder,pairedMSA_hhfilter_folder,pairedMSA_Nf90_folder):
    try:
        # here why use try???
        #here do not need to specifed correct path if in correct conda enviroment 
        hh_string="hhfilter -id 90 -cov 75 "+" -i " +pairedMSA_unfiltered_folder+\
        QueryPro_pro1+"and"+QueryPro_pro2+".fasta" +" -o "+\
        pairedMSA_hhfilter_folder+QueryPro_pro1+"and"+QueryPro_pro2+".A3M"+"  -M first"  # !!! here should be .A3M as it can only oouput this format file!! -M 50
        # although here we wrote it as fasta, but it actually A3M format 

        hh_string_escape=hh_string.replace("|","\|")
        #print(hh_string_escape)
        #print(hh_string)
        os.system(hh_string_escape)
    except:
        print("filtering problem ")
        print(QueryPro_pro1)
        print(QueryPro_pro2)

    hhfilter_file_name=pairedMSA_hhfilter_folder+QueryPro_pro1+"and"+QueryPro_pro2+".A3M"
    hhfilter_dict = SeqIO.to_dict(SeqIO.parse(hhfilter_file_name, "fasta"))
    hhfilter_keys=hhfilter_dict.keys()
    # why this step is necesseary 
    #Needed !!!! as after Nf90 threshold. we shoud use hhfilter paied MSA not original paired MSA
    #but hhfileter produced A2M format so we need to change it back to FSATA format !!!

    output_str=str()
    for name,record in paired_msa:
        if (name in hhfilter_keys):
            output_str=output_str+">"+name+"\n"+record+"\n"

    Nf90=len(hhfilter_keys)/(math.sqrt(QueryPro1_msa_len+QueryPro2_msa_len))
    
    if Nf90>Nf90_Thres:
        with open(pairedMSA_Nf90_folder+QueryPro_pro1+"and"+QueryPro_pro2+".fasta", "w") as myfile:
            myfile.write(output_str)
        
    
        return(Nf90)
    else:
        return(None)


def get_pairedMSA_inOneRun(row):
    currentSpe_TaxID,currentSpe_msa_removeGaps_path,QueryPro_pro1,QueryPro_pro2,Nf90_thres,pairedMSA_unfiltered_folder,pairedMSA_hhfilter_folder,pairedMSA_Nf90_folder=row

    # this check can be avoid now as we check this before this functin to avoid always trasmite large arument and check it many times 
    #if (QueryPro_pro1,QueryPro_pro2) not in pairedMSA_unfiltered_dict:


    # first find homologous form same species, 
    # then remove gaps. otherwise filtering on many sequenec will not used  at all
    # cancle this step as this will increase compuational time and i dont see much increase 
    getComStrSpeForPairedPro_Save_result=getComStrSpeForPairedPro_Save(currentSpe_TaxID,currentSpe_msa_removeGaps_path,QueryPro_pro1,QueryPro_pro2,pairedMSA_unfiltered_folder)
    #print("com_msa1 and com_msa2 len",len(common_msa1),len(common_msa2))
    if getComStrSpeForPairedPro_Save_result is not None: 
        QueryPro1_msa_len,QueryPro2_msa_len,fake_Nf90,paired_msa_list=getComStrSpeForPairedPro_Save_result


        if fake_Nf90>Nf90_thres:
            #print(fake_Nf90)
            Nf90=HHfilter2Save2Transfer2Fasta2Save(QueryPro_pro1,QueryPro_pro2,
                                                       QueryPro1_msa_len,QueryPro2_msa_len,
                                                       paired_msa_list,
                                                       Nf90_thres,
                                                       pairedMSA_unfiltered_folder,pairedMSA_hhfilter_folder,pairedMSA_Nf90_folder)
            #print(QueryPro_pro1,QueryPro_pro2,QueryPro1_msa_len,QueryPro2_msa_len,len(paired_msa_list),fake_Nf90,Nf90)
            if Nf90 is not None: 
                return([QueryPro_pro1,QueryPro_pro2,QueryPro1_msa_len,QueryPro2_msa_len,Nf90])
            else: 
                return (None)
        else:
            pass
            #print(fake_Nf90)

            
            
#heck paired MSA of homologus , if same proteins  from same species in both side
def getSameProteinRatio(record):
    pp1, pp2,final_msa_folder=record
    msa = AlignIO.read(final_msa_folder+pp1+"and"+pp2+".fasta", "fasta")

    sameProtein_count=0
    for record in msa:
        record_id=str(record.id)
        record_ids=record_id.split("and")
        if record_ids[0]==record_ids[1]:
            sameProtein_count +=1
    sameProtein_ratio=sameProtein_count/len(msa)
    return([pp1,pp2,sameProtein_ratio])
            
            

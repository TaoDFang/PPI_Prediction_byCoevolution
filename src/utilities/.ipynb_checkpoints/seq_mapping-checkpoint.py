import pickle
from Bio.Blast import NCBIXML

import subprocess
import os
import sys

def twoSepMapping(record):
    ''''''
    
    query_proID,subject_proID,overlap_method,query_fastaPath, subject_fastaPath,outputPath,blastp_path=record
    

    
    blastTableOutput_cmd = [
        blastp_path,
        '-query',query_fastaPath+query_proID+".fa",
        '-subject',subject_fastaPath+subject_proID+".fa",
        '-evalue',str(1e-6),# to make sure only meaningful hit are returned 
        '-outfmt', '6 qseqid  qlen sseqid  slen qstart qend sstart send evalue bitscore pident qcovs qcovhsp',
        '-out',outputPath+query_proID+"and"+subject_proID+".txt"


    ]



    skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.run(blastTableOutput_cmd, universal_newlines=True, **skw)
    p.returncode



    blastXMLOutput_cmd = [
        "blastp", 
        '-query',query_fastaPath+query_proID+".fa",
        '-subject',subject_fastaPath+subject_proID+".fa",
        '-evalue',str(1e-6),# to make sure only meaningful hit are returned 
        '-outfmt', str(5),
        '-out',outputPath+query_proID+"and"+subject_proID+".xml"


    ]



    skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.run(blastXMLOutput_cmd, universal_newlines=True, **skw)
    p.returncode

    with open( outputPath+query_proID+"and"+subject_proID+".xml") as result_handle:
        blast_record = NCBIXML.read(result_handle)



    subject_hsp_dict=dict()
    query_hsp_dict=dict()



    for alignment in blast_record.alignments:


        if overlap_method=="keep":
            sorted_alignment_hsps=sorted(alignment.hsps,key=lambda hsp: hsp.bits) 
        elif overlap_method=="remove":

            sorted_alignment_hsps=sorted(alignment.hsps,key=lambda hsp: hsp.bits,reverse=True) 

                
            best_hsp=sorted_alignment_hsps[0]
            keeped_subjectHSP_pos=[(best_hsp.sbjct_start,best_hsp.sbjct_end)] 
            keeped_queryHSP_pos=[(best_hsp.query_start,best_hsp.query_end)]  # notice here keeped_subjectHSP_pos and keeped_queryHSP_pos always have same length
            keeped_sorted_alignment_hsps=[best_hsp]
            for hsp in sorted_alignment_hsps[1:]: # check if include other hsps
                for idx in range(len(keeped_subjectHSP_pos)):
                    keeped_qstart,keeped_qend=keeped_subjectHSP_pos[idx]
                    keeped_sstart,keeped_send=keeped_queryHSP_pos[idx]

                    if ((hsp.sbjct_start>=keeped_qstart) and (hsp.sbjct_start<=keeped_qend)) or (hsp.sbjct_end>=keeped_qstart) and (hsp.sbjct_end<=keeped_qend):
                        pass
                    elif ((hsp.query_start>=keeped_sstart) and (hsp.query_start<=keeped_send)) or (hsp.query_end>=keeped_sstart) and (hsp.query_end<=keeped_send):
                        pass
                    else:
                        keeped_subjectHSP_pos.append((hsp.sbjct_start,hsp.sbjct_end))
                        keeped_queryHSP_pos.append((hsp.query_start,hsp.query_end))
                        keeped_sorted_alignment_hsps.append(hsp)
            # this step is not necesseary , just to keep consistent with "overlap_method=="keep""
            sorted_alignment_hsps=sorted(keeped_sorted_alignment_hsps,key=lambda hsp: hsp.bits) 

        for idx,hsp in enumerate(sorted_alignment_hsps):
                # notice here start and end position aer position in original , position of gaps are not included 
                subject_hsp_dict[idx]=(hsp.sbjct_start,hsp.sbjct_end,hsp.sbjct)
                query_hsp_dict[idx]=(hsp.query_start,hsp.query_end,hsp.query)




    # key idea here, if hsps overlap , then overlap it is 
    subject2query_AAmapDic=dict()
    for idx in range(len(subject_hsp_dict)):
        subject_hsp_tuple=subject_hsp_dict[idx]
        qstart,qend,qseq=subject_hsp_tuple

        query_hsp_tuple=query_hsp_dict[idx]
        sstart,send,sseq=query_hsp_tuple

        subject_gaps_soFar=0
        query_gaps_soFar=0
        for i, qlet in enumerate(qseq):
            slet=sseq[i] # for each hsp , when considering gaps, qseq and sseq always have same length 
            # and gaps from subject seq and sub seq cant be in same position
            # is it possible blast return two gaps in same position ??? even this is the case, my program can deal with this 
            if qlet != "-":
                # here only work if gaps is in subject and not in sub
                # this is true for current data is subject is also after trim ans thus shorter and thus need gaps
                #notice here postion are from balst and start from 1 not from 0
                if slet != "-":
                    subject2query_AAmapDic[qstart+i-subject_gaps_soFar]=sstart+i-query_gaps_soFar
                else: 
                    query_gaps_soFar +=1 # if qlet is not gap, but slet is gap,qlet residue cant be mapped 

            else:
                subject_gaps_soFar +=1  # if qlet is gap. then slet in this position cant be mapped 

                if slet =="-":
                    query_gaps_soFar +=1


    subject2query_AAmapDic={k:subject2query_AAmapDic[k] for k in sorted(subject2query_AAmapDic.keys())}
    #notice here postion are from balst and start from 1 not from 0
    with open(outputPath+query_proID+"and"+subject_proID+"subject2query_AAmapDic_overlapMethod"+overlap_method+".pickle", 'wb') as handle:
        pickle.dump(subject2query_AAmapDic, handle)


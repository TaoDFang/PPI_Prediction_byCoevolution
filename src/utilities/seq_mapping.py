import pickle
from Bio.Blast import NCBIXML

import subprocess

# code adapted  from /code/MNF/src/tao_utilities/protein_otherFeatures.py
# and /code/MNF/notebooks/Ecoli_PDB/Physical_Distance_Calculation.ipynb


#here chante the name query and subject in old script to the subject and query (reverse it to fit to  the usually customer)
def twoSepMapping(record):
    '''this function was also used in /code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_intergrationAtResidueLevel.ipynb'''
    
    query_proID,subject_proID,overlap_method,query_fastaPath, subject_fastaPath,outputPath,blastp_path=record
    

    
    blastTableOutput_cmd = [
        blastp_path,
        '-query',query_fastaPath+query_proID+".fa",
        '-subject',subject_fastaPath+subject_proID+".fa",
        '-evalue',str(1e-6),# to make sure only meaningful hit are returned 
        '-outfmt', '6 qseqid  qlen sseqid  slen qstart qend sstart send evalue bitscore pident qcovs qcovhsp',
        '-out',outputPath+query_proID+"and"+subject_proID+".txt"


    ]
    
    # print("    ".join(blastTableOutput_cmd))


    skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.run(blastTableOutput_cmd, universal_newlines=True, **skw)
    p.returncode



    blastXMLOutput_cmd = [
        "/mnt/mnemo5/tao/BeeiveProgram/ncbi-blast-2.10.0+/bin/blastp", 
        '-query',query_fastaPath+query_proID+".fa",
        '-subject',subject_fastaPath+subject_proID+".fa",
        '-evalue',str(1e-6),# to make sure only meaningful hit are returned 
        '-outfmt', str(5),
        '-out',outputPath+query_proID+"and"+subject_proID+".xml"


    ]



    skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.run(blastXMLOutput_cmd, universal_newlines=True, **skw)
    p.returncode


    # use totorial from biopython https://biopython.readthedocs.io/en/latest/chapter_blast.html
    with open( outputPath+query_proID+"and"+subject_proID+".xml") as result_handle:
        blast_record = NCBIXML.read(result_handle)


    # one problem here is that if we process by hsp. probabily not all resiiue in queyr protein could be mapped
    # which is okay,as when two protein are not have same len, of course some residues in short protein  could not be mapped 
    subject_hsp_dict=dict()
    query_hsp_dict=dict()

    # ???? here there is a problem is that if blast_record.alignments is None
    # then this for loop will never get started 
    # lucky here is oaky as we want subject_hsp_dict and query_hsp_dict are empty dictionary in this case 

    for alignment in blast_record.alignments:
#         if len(alignment.hsps)>2:
#             print(query_proID,subject_proID,query_fastaPath, subject_fastaPath,outputPath)

        if overlap_method=="keep":
            #here sorted hsps by thier bits scores from smallest to largest
            # so if these hsps overlaps , the residue mapping realtion in best hsp are used 
            sorted_alignment_hsps=sorted(alignment.hsps,key=lambda hsp: hsp.bits) 
        elif overlap_method=="remove":
            # there sorted hsp by their bits scores from largest to smallest
            # if there are overlap. remove hsp with lower bits score
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

#         if (overlap_method=="remove") and (len(sorted_alignment_hsps)>1):
#             print(query_proID,subject_proID,query_fastaPath, subject_fastaPath,outputPath)

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
                # this is true for current data is subject is alos after trim ans thus shorter and thus need gaps

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
    #return(subject2query_AAmapDic)


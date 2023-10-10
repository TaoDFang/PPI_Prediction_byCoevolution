import os 
import pandas as pd 
import numpy as np 
import subprocess
import re 
from Bio import AlignIO 
import pickle

def ParseCurSpeFastaByProteins_STRNG1105(currentSpe_fastaData,currentSpeProSeqPath_ByProteins):
    print("currentSpe_fastaData,currentSpeProSeqPath_ByProteins:",currentSpe_fastaData,currentSpeProSeqPath_ByProteins)
    with open(currentSpe_fastaData, "r") as inFile:
        for line in inFile:
            line = line.rstrip() 
            if line.startswith('>'):
                pro_name=line[1:]
                with open(currentSpeProSeqPath_ByProteins+pro_name + ".fa", 'a') as outFile:
                     outFile.write('%s\n' % line)

            else:
                with open(currentSpeProSeqPath_ByProteins+pro_name + ".fa", 'a') as outFile:
                    outFile.write("%s\n" % line)
                    
                    
def combinedRBH_STRING1105(record):
    file, allBacteriaProteomeNum=record
    base_file=os.path.basename(file)
    base_file=base_file[:-10]

    currentSpe_ogs=pd.read_csv(file,
                                         header=None,index_col=None,sep="\t",dtype={0: str,1:str})
    #if currentSpe_ogs.shape[0]>400: # here 400 was choose as 1% of all proteome used in orignal david bakers paper
    if currentSpe_ogs.shape[0]>(allBacteriaProteomeNum*0.01):
        
        currentSpe_ogs_pros=currentSpe_ogs.iloc[1:,1].values.tolist()
    

        overlap_pros=currentSpe_ogs_pros
        new_frame=pd.DataFrame({0:[base_file]*len(overlap_pros),1:list(overlap_pros)})

        return(new_frame)
    
    
    
def fun_newSingleMSA_EggNOG_OrthologousGroup_faidx(record):
    # first need to make sure script is executable : chmod +x name_of_your_file_script
    code_utilities_folder,currentSpe_TaxID,OG_idx,currentSpe_fastaData,origSTRINGBacteriaProSeqPath,currentSpe_OrthologousGroup_Fa_path,currentSpe_OrthologousGroup_Fa_logpath,newsingleMSA_RBH_OrthologousGroup_fileName=record
 
    
    log_name = os.path.join(currentSpe_OrthologousGroup_Fa_logpath, str(OG_idx) + '.log')
    with open(log_name, 'w') as logfp:
        cmd = ["%sbash_scripts/newSingleMSA_EggNOG_OrthologousGroup_faidx.sh %s %s  %s %s %s %s" % (code_utilities_folder,currentSpe_TaxID,str(OG_idx),currentSpe_fastaData,origSTRINGBacteriaProSeqPath,currentSpe_OrthologousGroup_Fa_path,newsingleMSA_RBH_OrthologousGroup_fileName),
              ] 
        skw = dict(stdout=logfp, stderr=logfp)

        p = subprocess.run(cmd, universal_newlines=True,shell=True, **skw,) # add shell=True to run bash script 
        p.check_returncode()
    
    
def fun_newSingleMSA_EggNOG_OrthologousGroup_phmmer(record):
    # first need to make sure script is executable : chmod +x name_of_your_file_script
    code_utilities_folder,phmmer_path,idx,currentSpe_OrthologousGroup_Fa_path,currentSpeProSeqPath_ByProteins,currentSpe_phmmer_outPath,currentSpe_phmmer_logPath=record
    
    log_name = os.path.join(currentSpe_phmmer_logPath, str(idx) + '.log')
    
    with open(log_name, 'w') as logfp:
        
        cmd = ["%sbash_scripts/newSingleMSA_EggNOG_OrthologousGroup_phmmer.sh %s %s  %s %s %s" % (code_utilities_folder,phmmer_path,str(idx),currentSpe_OrthologousGroup_Fa_path,currentSpeProSeqPath_ByProteins,currentSpe_phmmer_outPath),

              ] 
        skw = dict(stdout=logfp, stderr=logfp)

        p = subprocess.run(cmd, universal_newlines=True,shell=True, **skw,) # add shell=True to run bash script 
        p.check_returncode()
    
    
    
def ClustoMSA(record):
    name,currentSpe_phmmer_OrthologousGroup_path,currentSpe_ClustoMSA_path,clustalo_path=record
    cmd = [
        clustalo_path,
        '-i', currentSpe_phmmer_OrthologousGroup_path+name+".fa",
        '-o', currentSpe_ClustoMSA_path+name,
        "--outfmt=st","--threads=1"

    ]

    skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.run(cmd, universal_newlines=True, **skw)
    p.check_returncode()
    
    
    
def hmmbuild(record):
    name,currentSpe_hmm_profiles_path,currentSpe_hmm_profiles_path,currentSpe_ClustoMSA_path,hmmbuild_path=record
    cmd = [
        hmmbuild_path, 
        '--cpu', "1",
        '-o', currentSpe_hmm_profiles_path+name+".outinfo",
        currentSpe_hmm_profiles_path+name+".hmm",
        currentSpe_ClustoMSA_path+name

    ]
    skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.run(cmd, universal_newlines=True, **skw)
    p.check_returncode()
    
    
def hmmalign(record):
    name,currentSpe_hmm_align_path,currentSpe_hmm_profiles_path,currentSpe_OrthologousGroup_Fa_path,hmmalign_path=record
    cmd = [hmmalign_path, 
        '-o', currentSpe_hmm_align_path+name+".sto",
        currentSpe_hmm_profiles_path+name+".hmm",
        currentSpe_OrthologousGroup_Fa_path+name+".fa"

    ]

    skw = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p = subprocess.run(cmd, universal_newlines=True, **skw)
    p.check_returncode()
    

    
    
#here track position of aa of seq in MSA in original protein seq 
def removeGapsANDtrackAAPosOfSeqInMSA(record):
    #file_name=hmmalign_sto_files[idx]
    #print(file_name)
    file_name,currentSpe_ClustoMSA_path,currentSpe_msa_trackGapsPos_path,currentSpe_msa_removeGaps_path=record
    gap_seq_ratio=0.5,
    gap_col_ratio=0.5
    ################
    alphabet = "ARNDCQEGHILKMFPSTWYV-"
    states = len(alphabet)

    a2n = {}
    for a,n in zip(alphabet,range(states)):
        a2n[a] = n
    def aa2num(aa):
        if aa in a2n: 
            return a2n[aa]
        else: 
            #print(aa)
            return a2n['-']

    n2a = {}
    for k , v in a2n.items():
        n2a[v]=k

    reTemp=re.compile(".*\/(.*)\.sto")
    currentSpe_pro_name=reTemp.match(file_name).group(1)
    #print(currentSpe_pro_name)
    sto_align = AlignIO.read(file_name, "stockholm")

    header = []
    sequence_digit = []

    for sto_record in sto_align:
        sto_id=sto_record.id
        sto_seq=str(sto_record.seq)
        sto_seq=str(sto_seq).upper()
        sto_seq=sto_seq.replace(".","-")
        header.append(sto_id)
        sequence_digit.append([aa2num(aa) for aa in sto_seq])

    header_array=np.array(header)
    seq_array=np.array(sequence_digit)


    # deal with ref seq, and remove gap in ref
    ref_tmp = seq_array[0,:]
    ref_gaps = np.where(ref_tmp==(states-1))[0]
    #print("len(ref_gaps):",len(ref_gaps))
    ref_non_gaps = np.where(ref_tmp!=(states-1))[0]



    # remove gap in seq 
    # this step affect positive tracking 
    ref_seq_array=seq_array[:,ref_non_gaps]
    # here is not wrong even if we incldue ref-seq, as gaps in ref-seq have been removed in last step 
    seq_tmp = (ref_seq_array == states-1).astype(np.float)

    # this  could make jupyter notebook dead 
    
    if ref_seq_array.shape[1]!=0:
        seq_non_gaps = np.where(np.sum(seq_tmp,1)/ref_seq_array.shape[1] < gap_seq_ratio)[0]
    else: 
        Clusto_msa=AlignIO.read(currentSpe_ClustoMSA_path+currentSpe_pro_name,"stockholm")
        print("Ref seq are all gaps",file_name,seq_array.shape,ref_seq_array.shape)
        print("Seed alignment: ",currentSpe_ClustoMSA_path+currentSpe_pro_name, Clusto_msa.get_alignment_length())
        return 


    # remove gap in column
    # here is important it use seq_array again  insted of ref_seq_array
    # as column which  is gap  in ref-seq are kept
    # so that  position  of gaps is the position in origial seq which 
    # is  esier position  tracking 
    # as here gaps in column will also include some gaps in ref-seq 
    


    seq_seq_array=seq_array[seq_non_gaps,:]
    seq_header_array=header_array[seq_non_gaps]
    #print(seq_array.shape,seq_seq_array.shape)
    tmp = (seq_seq_array == states-1).astype(np.float)
    # its fine to fine to keep index here to late comparing reasone as same gap position in ref seq will be removed 
    # here this will not affect non gaps in first row 
    # and for gaps in first row, it will be removed before anyway 
    if seq_seq_array.shape[0]!=0:
        col_gaps = np.where(np.sum(tmp,0)/seq_seq_array.shape[0] > gap_col_ratio)[0]
    else:
        print("All seqs have too many gaps",file_name,seq_array.shape,seq_seq_array.shape)
        return 
        
    #print("len(col_gaps):",len(col_gaps))
    col_gaps=[g for g in col_gaps if g  not in ref_gaps]
    #print("len(col_gaps):",len(col_gaps))


    ref_col_gaps=sorted(list(set(ref_gaps.tolist()+col_gaps)))
    #print("len(ref_col_gaps)",len(ref_col_gaps))



    # get positive of aa of proteins in original seq 
    # check aa is not gap and not remove in gap filtering step in ref-seq and in all colunns

    allSeq_keptAA_posInOrigalPro=list()
    for i in range(seq_header_array.shape[0]):
        seq_id=seq_header_array[i]
        seq=seq_seq_array[i,:] # here seq includ gaps, not examly original protein seq 

        keptAA_posInOrigalPro=list()
        original_pos=0
        for i,s in enumerate(seq): 

            if (s !=(states-1)) and (i not in ref_gaps) and (i not in col_gaps):
                keptAA_posInOrigalPro.append(original_pos)

            if (s !=(states-1)):
                original_pos +=1

        allSeq_keptAA_posInOrigalPro.append([seq_id,keptAA_posInOrigalPro])
        # here keptAA_posInOrigalPro is a list containing position of each AA(of final MSA ) in each  origin seq 
        # here seq_seq_array is MSA contain gaps in columns 
        # remember for final MSA after all gaps removeing , first row contain no gaps  
        
        
        # key/index is position in orig seq ,v/value is pos in msa seq  (with gaps removed )
        # consider we count original_pos as long as its not gap 
        
        # no no no no no key/index is position in orig seq ,v/value is pos in msa seq  (with gaps removed )
        # keep in mind, pos in orig seq  is shorter
        # remember for paired MSA, first row/currentSpe_pro seq, all gaps shou be already removed 
        # but for others rows, we need to consider gaps 
        
    with open(currentSpe_msa_trackGapsPos_path+currentSpe_pro_name+'_allSeq_keptAA_posInOrigalPro.pickle', 'wb') as handle:
        pickle.dump(allSeq_keptAA_posInOrigalPro, handle)
    
    
    # get and save final MSA  to fasta format 

    final_seq_array=np.delete(seq_seq_array,ref_col_gaps,1)
    final_header_array=seq_header_array
    
    sto_outStr=""      
    for i in range(final_header_array.shape[0]):
        seq_digit=final_seq_array[i,:]
        seq=[n2a[n] for n in seq_digit]
        sto_outStr=sto_outStr+">"+final_header_array[i]+"\n"+"".join(seq)+"\n"

    with open(currentSpe_msa_removeGaps_path+currentSpe_pro_name+".fa","w") as outfile:
        outfile.write(sto_outStr)   
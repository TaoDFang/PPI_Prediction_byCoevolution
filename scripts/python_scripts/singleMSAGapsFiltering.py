import argparse
import os 
import glob
import multiprocessing as mp
from Bio import AlignIO 
import pickle
import math

from create_singleMSA import removeGapsANDtrackAAPosOfSeqInMSA

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--currentSpe_hmm_align_path', type=str, help='currentSpe_hmm_align_path')
    parser.add_argument('-c','--currentSpe_ClustoMSA_path', type=str, help='currentSpe_ClustoMSA_path')
    parser.add_argument('-tp','--currentSpe_msa_trackGapsPos_path', type=str, help='currentSpe_msa_trackGapsPos_path')
    parser.add_argument('-rp','--currentSpe_msa_removeGaps_path', type=str, help='currentSpe_msa_removeGaps_path')
    parser.add_argument('-m','--currentSpeMSAGapsFilteringMetaFolder', type=str, help='currentSpeMSAGapsFilteringMetaFolder')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')


    args = parser.parse_args()
    currentSpe_hmm_align_path=args.currentSpe_hmm_align_path
    currentSpe_ClustoMSA_path=args.currentSpe_ClustoMSA_path
    currentSpe_msa_trackGapsPos_path=args.currentSpe_msa_trackGapsPos_path
    currentSpe_msa_removeGaps_path=args.currentSpe_msa_removeGaps_path
    currentSpeMSAGapsFilteringMetaFolder=args.currentSpeMSAGapsFilteringMetaFolder
    mp_task_nums=int(args.mp_task_nums)
    
    
    


    hmmalign_sto_files=glob.glob(os.path.join(currentSpe_hmm_align_path,"*.sto"))
    print("len(hmmalign_sto_files):",len(hmmalign_sto_files))
    hmmalign_sto_files_ArgForremoeGaps=[(file_name,currentSpe_ClustoMSA_path,currentSpe_msa_trackGapsPos_path,currentSpe_msa_removeGaps_path) for file_name in hmmalign_sto_files]
    pool = mp.Pool(mp_task_nums)
    pool.map(removeGapsANDtrackAAPosOfSeqInMSA, hmmalign_sto_files_ArgForremoeGaps)
    pool.close()
    # CPU times: user 3.73 s, sys: 3.18 s, total: 6.9 s
    # Wall time: 17min 22s


    # # get final single MSA len
    final_MSAs=glob.glob(currentSpe_msa_removeGaps_path+"*")
    print("len(final_MSAs):",len(final_MSAs))
    print(final_MSAs[0])

    final_MSA_pros=[os.path.basename(f) for f in final_MSAs]
    final_MSA_pros=[f[:-3] for f in final_MSA_pros]
    print(final_MSA_pros[0])

    # get len and Nf90 of all single MSA 

    fasta_protein_lens=dict()
    fasta_protein_Nf90s=dict()

    for pro in final_MSA_pros:
        alig = AlignIO.read(currentSpe_msa_removeGaps_path+pro+".fa", "fasta")
        fasta_protein_lens[pro]=len(alig[0])
        fasta_protein_Nf90s[pro]=len(alig)/math.sqrt(len(alig[0]))

    with open(currentSpeMSAGapsFilteringMetaFolder+'fasta_protein_lens_dict.pickle', 'wb') as handle:
        pickle.dump(fasta_protein_lens, handle)

    with open(currentSpeMSAGapsFilteringMetaFolder+'fasta_protein_Nf90s_dict.pickle', 'wb') as handle:
        pickle.dump(fasta_protein_Nf90s, handle)



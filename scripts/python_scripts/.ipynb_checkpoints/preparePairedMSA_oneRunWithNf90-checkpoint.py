import argparse
import os 
import glob
import multiprocessing as mp
import pandas as pd 
from multiprocessing import get_context
import pickle

from  create_pairedMSA import get_pairedMSA_inOneRun


import random 
random.seed(10)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--currentSpe_TaxID', type=str, help='currentSpe_TaxID')
    parser.add_argument('-m','--currentSpeMiddleDataPath', type=str, help='currentSpeMiddleDataPath')
    parser.add_argument('-rp','--currentSpe_msa_removeGaps_path', type=str, help='currentSpe_msa_removeGaps_path') # notice currentSpe_msa_removeGaps_path is the old currentSpe_hmmalign_path in old script ,here give a differnt name for consitency of previous variable names 
    parser.add_argument('-un','--pairedMSA_unfiltered_folder', type=str, help='pairedMSA_unfiltered_folder')
    parser.add_argument('-hh','--pairedMSA_hhfilter_folder', type=str, help='pairedMSA_hhfilter_folder')
    parser.add_argument('-nf90f','--pairedMSA_Nf90_folder', type=str, help='pairedMSA_Nf90_folder')
    parser.add_argument('-nf90csv','--pairedMSA_Nf90_csv', type=str, help='pairedMSA_Nf90_csv')
    parser.add_argument('-nf90','--Nf90_thres', type=str, help='Nf90_thres')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')

    args = parser.parse_args()
    currentSpe_TaxID=args.currentSpe_TaxID
    currentSpeMiddleDataPath=args.currentSpeMiddleDataPath
    currentSpe_msa_removeGaps_path=args.currentSpe_msa_removeGaps_path
    pairedMSA_unfiltered_folder=args.pairedMSA_unfiltered_folder
    pairedMSA_hhfilter_folder=args.pairedMSA_hhfilter_folder
    pairedMSA_Nf90_folder=args.pairedMSA_Nf90_folder
    pairedMSA_Nf90_csv=args.pairedMSA_Nf90_csv
    Nf90_thres=int(args.Nf90_thres)
    mp_task_nums=int(args.mp_task_nums)


    with open(currentSpeMiddleDataPath+'fasta_protein_lens_dict.pickle', 'rb') as handle:
        fasta_protein_lens=pickle.load(handle)

    with open(currentSpeMiddleDataPath+'fasta_protein_Nf90s_dict.pickle', 'rb') as handle:
        fasta_protein_Nf90s=pickle.load(handle)


    currentSpe_final_MSAs=glob.glob(currentSpe_msa_removeGaps_path+"*")
    print("len(currentSpe_final_MSAs):",len(currentSpe_final_MSAs))
    print(currentSpe_final_MSAs[0])

    currentSpe_final_pros=[os.path.basename(f) for f in currentSpe_final_MSAs]
    currentSpe_final_pros=[p[:-3] for p in currentSpe_final_pros]
    print(currentSpe_final_pros[0])

    currentSpe_final_pro_dict=dict([(p,1)for p in currentSpe_final_pros])



    if not os.path.exists(pairedMSA_Nf90_csv):
        with open(pairedMSA_Nf90_csv, 'w') as fp:
            pass


    fakeNf90_fasta_proteins=[k for k , v in fasta_protein_Nf90s.items() if v>=Nf90_thres]
    fakeNf90_fasta_proteins=sorted(fakeNf90_fasta_proteins)
    print("len(fakeNf90_fasta_proteins):",len(fakeNf90_fasta_proteins))
    print(fakeNf90_fasta_proteins[0:3])

    print("len(fakeNf90_fasta_proteins):",len(fakeNf90_fasta_proteins))

    allPPIs=[(fakeNf90_fasta_proteins[i],fakeNf90_fasta_proteins[j]) for i in range(len(fakeNf90_fasta_proteins)-1) for j in range(i+1,len(fakeNf90_fasta_proteins)) ]
    print("len(allPPIs):",len(allPPIs),len(fakeNf90_fasta_proteins)*(len(fakeNf90_fasta_proteins)-1)/2)
    print(allPPIs[0:3])

    random.seed(10)
    random.shuffle(allPPIs)
    print("len(allPPIs):",len(allPPIs),len(fakeNf90_fasta_proteins)*(len(fakeNf90_fasta_proteins)-1)/2)
    print(allPPIs[0:3])

    #%%time 
    # now = datetime.now()
    # print("detect existed pairedMSA_unfiltered_files begins: ", now)

    pairedMSA_unfiltered_files =os.listdir(pairedMSA_unfiltered_folder) # here to check if paired MSA already been processed before , no use in nextflow as everytime the process start from null output 
    print("len(pairedMSA_unfiltered_files):",len(pairedMSA_unfiltered_files))
    print(pairedMSA_unfiltered_files[0:3])

    pairedMSA_unfiltered_pps=[f.split("and") for f in pairedMSA_unfiltered_files]
    pairedMSA_unfiltered_pps=[(p1,p2[0:-6]) for p1,p2 in pairedMSA_unfiltered_pps]

    pairedMSA_unfiltered_dict=dict([(f,1) for f in pairedMSA_unfiltered_pps])


    # now = datetime.now()
    # print("detect existed pairedMSA_unfiltered_files ends: ", now) 
    # this step takes around 6 minutes when run in notebook in phobos
    # while it takes around 1 hours when run with script in phobos 

    #%%time
    currentSpe_allPPIs_forPairMSA_dict={pp:1 for pp in allPPIs if pp not in pairedMSA_unfiltered_dict}
    currentSpe_allPPIs_ArgForPairMSA=[(currentSpe_TaxID,currentSpe_msa_removeGaps_path,p1,p2,Nf90_thres,pairedMSA_unfiltered_folder,pairedMSA_hhfilter_folder,pairedMSA_Nf90_folder) for p1,p2 in currentSpe_allPPIs_forPairMSA_dict]
    #this step for debug reason, comment it after make sure no bugs
    currentSpe_allPPIs_ArgForPairMSA=currentSpe_allPPIs_ArgForPairMSA[0:100]


    print("len(currentSpe_allPPIs_ArgForPairMSA):",len(currentSpe_allPPIs_ArgForPairMSA))



    if len(currentSpe_allPPIs_ArgForPairMSA)>0:
        #pool=mp.Pool(50)   #mp.Pool(160) 
        pool=get_context("spawn").Pool(mp_task_nums)
        get_pairedMSA_inOneRun_result=pool.map(get_pairedMSA_inOneRun,currentSpe_allPPIs_ArgForPairMSA)
        pool.close() 

        get_pairedMSA_inOneRun_result=[s for s in get_pairedMSA_inOneRun_result if s is not None]
        get_pairedMSA_inOneRun_result_frame=pd.DataFrame(get_pairedMSA_inOneRun_result,columns=["protein1","protein2","L1","L2","Nf90"])
        #in the nextflow piple, the mode ="a" is meaningless as everytime the process start from begining and preprall all possible ppi 
        get_pairedMSA_inOneRun_result_frame.to_csv(pairedMSA_Nf90_csv,mode="a",
                            header=None,index=None,sep="\t")
        print("get_pairedMSA_inOneRun_result_frame.shape:",get_pairedMSA_inOneRun_result_frame.shape)

  

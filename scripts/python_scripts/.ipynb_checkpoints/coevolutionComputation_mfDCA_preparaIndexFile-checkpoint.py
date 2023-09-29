import argparse
import os 
import glob
import csv
import multiprocessing as mp
from multiprocessing import get_context
import pickle
from datetime import datetime




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-dpath','--DCA_coevolutoin_path', type=str, help='DCA_coevolutoin_path')
    parser.add_argument('-dipath','--IndexDCA_coevolutoin_path', type=str, help='IndexDCA_coevolutoin_path')
    parser.add_argument('-m','--currentSpeMSAGapsFilteringMetaFolder', type=str, help='currentSpeMSAGapsFilteringMetaFolder')
    parser.add_argument('-acsv','--PPIInfoBeforeCoEvoComp_csv', type=str, help='PPIInfoBeforeCoEvoComp_csv')
    parser.add_argument('-nf90f','--pairedMSA_Nf90_folder', type=str, help='pairedMSA_Nf90_folder')
    parser.add_argument('-bs','--blockSize', type=str, help='blockSize')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    
    
    args = parser.parse_args()
    currentSpeMSAGapsFilteringMetaFolder=args.currentSpeMSAGapsFilteringMetaFolder
    PPIInfoBeforeCoEvoComp_csv=args.PPIInfoBeforeCoEvoComp_csv
    DCA_coevolutoin_path=args.DCA_coevolutoin_path
    IndexDCA_coevolutoin_path=args.IndexDCA_coevolutoin_path
    pairedMSA_Nf90_folder=args.pairedMSA_Nf90_folder
    blockSize=int(args.blockSize)
    mp_task_nums=int(args.mp_task_nums)

    
    with open(currentSpeMSAGapsFilteringMetaFolder+'fasta_protein_lens_dict.pickle', 'rb') as handle:
        fasta_protein_lens=pickle.load(handle)

    # use  postive and negative ppi to only use pp that have large Nf90 value, and after removing deep homologs  
    # currentSpe_allPPIs_beforeCoEvoComp_frame = pd.read_csv(PPIInfoBeforeCoEvoComp_csv,header=0,index_col=None,sep="\t")
    # print("currentSpe_allPPIs_beforeCoEvoComp_frame.shape:",currentSpe_allPPIs_beforeCoEvoComp_frame.shape)
    # currentSpe_allPPIs_beforeCoEvoComp_info=currentSpe_allPPIs_beforeCoEvoComp_frame.values.tolist()
    currentSpe_allPPIs_beforeCoEvoComp_info=list()
    with open(PPIInfoBeforeCoEvoComp_csv, "r") as file:
        my_reader = csv.reader(file, delimiter="\t")
        next(my_reader, None)  # skip the headers
        for row in my_reader:
            currentSpe_allPPIs_beforeCoEvoComp_info.append(row)

            
    currentSpe_allPPIs_pps=[(p1,p2) for p1, p2,_ in currentSpe_allPPIs_beforeCoEvoComp_info]
    print("len(currentSpe_allPPIs_pps):",len(currentSpe_allPPIs_pps))
    


    #if re-run, first check if file is already existed, need to change pythion code 
    existed_pydcaFNAPC_files=glob.glob(DCA_coevolutoin_path+"*_pydcaFNAPC_array.npz")
    print("len(existed_pydcaFNAPC_files):",len(existed_pydcaFNAPC_files))
    print("existed_pydcaFNAPC_files[0:3]:",existed_pydcaFNAPC_files[0:3])
    existed_pydcaFNAPC_files=[os.path.basename(f) for f in existed_pydcaFNAPC_files]
    existed_pydcaFNAPC_pps=[f.split("and") for f in existed_pydcaFNAPC_files]
    existed_pydcaFNAPC_pps=[(p1,p2[:-21]) for p1, p2 in existed_pydcaFNAPC_pps]
    existed_pydcaFNAPC_pp_dict=dict([(pp,1) for pp in existed_pydcaFNAPC_pps])



    #%%time 
    currentSpe_allPPIs_pps_forDCA=[pp for pp in currentSpe_allPPIs_pps if pp not in existed_pydcaFNAPC_pp_dict]
    currentSpe_allPPIs_pps_ArgForDCA=[(p1,p2,fasta_protein_lens[p1],fasta_protein_lens[p2],pairedMSA_Nf90_folder,DCA_coevolutoin_path) for p1 , p2 in currentSpe_allPPIs_pps_forDCA]
    #for dedub reson, delete it later 
    currentSpe_allPPIs_pps_ArgForDCA=currentSpe_allPPIs_pps_ArgForDCA[0:3]
    print("len(currentSpe_allPPIs_pps_ArgForDCA):",len(currentSpe_allPPIs_pps_ArgForDCA))



    # # here no need anymore, in nextlfow process , each time IndexDCA_coevolutoin_path is a contempory folder 
    # IndexDCA_coevolutoin_path_files = glob.glob(IndexDCA_coevolutoin_path+"*.csv")
    # for f in IndexDCA_coevolutoin_path_files:
    #     os.remove(f)

    # blockSize=500
    index_count=0
    for bidx in range(0,len(currentSpe_allPPIs_pps_forDCA),blockSize):
        block_currentSpe_allPPIs_pps_forDCA=currentSpe_allPPIs_pps_forDCA[bidx:(bidx+blockSize)]
        #print(len(block_currentSpe_allPPIs_pps_forDCA))
        block_currentSpe_allPPIs_pps_forDCA_frame=pd.DataFrame(block_currentSpe_allPPIs_pps_forDCA)
        block_currentSpe_allPPIs_pps_forDCA_frame.to_csv(IndexDCA_coevolutoin_path+str(index_count)+".csv",
                                                                header=None,index=None,sep="\t")
        index_count=index_count+1


    print("index_count:",index_count)


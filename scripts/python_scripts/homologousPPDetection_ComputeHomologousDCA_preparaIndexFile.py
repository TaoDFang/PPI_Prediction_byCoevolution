import argparse
import os 
import glob
import pandas as pd 
import multiprocessing as mp
from multiprocessing import get_context
import pickle
from datetime import datetime




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-q','--Query_tuple', type=str, help='Query_tuple')
    parser.add_argument('-s','--Subject_tupleList', type=str, help='Subject_tupleList')
    parser.add_argument('-f','--newSTRING_rootFolder', type=str, help='newSTRING_rootFolder')
    parser.add_argument('-c','--CoEvo_data_folder', type=str, help='CoEvo_data_folder')
    parser.add_argument('-bn','--blockNum', type=str, help='blockNum')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    
    args = parser.parse_args()
    Query_tuple=args.Query_tuple
    Subject_tupleList=args.Subject_tupleList
    newSTRING_rootFolder=args.newSTRING_rootFolder
    CoEvo_data_folder=args.CoEvo_data_folder
    blockNum=int(args.blockNum)
    mp_task_nums=int(args.mp_task_nums)
    
    Query_tuple=Query_tuple.split("_")
    
    Subject_tupleList=Subject_tupleList.split("_")
    Subject_tupleList=[(Subject_tupleList[i],Subject_tupleList[i+1]) for i in range(0,len(Subject_tupleList)-1,2)]
    
    
    for currentSubject_EggNOG_maxLevel,currentSubject_TaxID in Subject_tupleList:
        print(currentSubject_EggNOG_maxLevel,currentSubject_TaxID)
        
        Query_prefix="BestHomologousPPFor"+Query_tuple[1]+"AtEggNOGmaxLevel"+Query_tuple[0]+"_"
        
        input_root_folder=CoEvo_data_folder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_eggNOGfilteredData/"
        pairedMSA_Nf90_folder = input_root_folder+"pair_MSA_Nf90_PasteAlign/"
        DCA_coevolutoin_path=input_root_folder+"coevolutoin_result_DCA/"
        IndexDCA_coevolutoin_path=input_root_folder+"coevolutoin_computation_IndexDCA/"
        if not os.path.exists(DCA_coevolutoin_path):
            os.makedirs(DCA_coevolutoin_path)
            
        if not os.path.exists(IndexDCA_coevolutoin_path):
            os.makedirs(IndexDCA_coevolutoin_path)
            
            
        currentSpeMSAGapsFilteringMetaFolder=newSTRING_rootFolder+currentSubject_TaxID+"_EggNOGmaxLevel"+currentSubject_EggNOG_maxLevel+"_MSAGapsFiltering/"
        with open(currentSpeMSAGapsFilteringMetaFolder+'fasta_protein_lens_dict.pickle', 'rb') as handle:
            fasta_protein_lens=pickle.load(handle)

        
        Benchmark_folder=input_root_folder+Query_prefix+"AllPPI_Benchmark/" #
        Current_Subject_homologousPPs_beforeDCAcomputation_file=Benchmark_folder+"Current_Subject_BestHomologousPPs_beforeDCAcomputation.pickle"
        with open(Current_Subject_homologousPPs_beforeDCAcomputation_file, 'rb') as handle:
            Current_Subject_homologousPPs=pickle.load(handle)


        # ### DCA


        #if re-run, first check if file is already existed, need to change pythion code 

        existed_pydcaFNAPC_files=glob.glob(DCA_coevolutoin_path+"*_pydcaFNAPC_array.npz")
        print(len(existed_pydcaFNAPC_files))
        print(existed_pydcaFNAPC_files[0:3])

        existed_pydcaFNAPC_files=[os.path.basename(f) for f in existed_pydcaFNAPC_files]
        existed_pydcaFNAPC_pps=[f.split("and") for f in existed_pydcaFNAPC_files]
        existed_pydcaFNAPC_pps=[(p1,p2[:-21]) for p1, p2 in existed_pydcaFNAPC_pps]

        existed_pydcaFNAPC_pp_dict=dict([(pp,1) for pp in existed_pydcaFNAPC_pps])



        #%%time 
        currentSubject_allPPI_forDCA=[pp for pp in Current_Subject_homologousPPs if pp not in existed_pydcaFNAPC_pp_dict] 
        #here 0:3 for testing reson , comment it laterß
        currentSubject_allPPI_ArgForDCA=[(p1,p2,fasta_protein_lens[p1],fasta_protein_lens[p2],pairedMSA_Nf90_folder,DCA_coevolutoin_path) for p1 , p2 in currentSubject_allPPI_forDCA[0:6]]
        print(len(currentSubject_allPPI_ArgForDCA))



        # unlike for query species, for subject species, here is necesseary 
        IndexDCA_coevolutoin_path_files = glob.glob(IndexDCA_coevolutoin_path+"*.csv")
        for f in IndexDCA_coevolutoin_path_files:
            os.remove(f)

        assert blockNum<len(currentSubject_allPPI_ArgForDCA)

        block_range=np.linspace(0,len(currentSubject_allPPI_ArgForDCA),blockNum+1,endpoint=True,dtype=int)
        for i in range(blockNum):
            block_currentSpe_allPPIs_pps_forDCA=currentSubject_allPPI_ArgForDCA[block_range[i]:block_range[i+1]]
            block_currentSpe_allPPIs_pps_forDCA_frame=pd.DataFrame(block_currentSpe_allPPIs_pps_forDCA)
            block_currentSpe_allPPIs_pps_forDCA_frame.to_csv(IndexDCA_coevolutoin_path+str(i)+".csv",
                                                                    header=None,index=None,sep="\t")  
        
import argparse
import os 
import glob
import multiprocessing as mp

from create_singleMSA import fun_newSingleMSA_EggNOG_OrthologousGroup_phmmer


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-b','--currentSpeProSeqPath_ByProteins', type=str, help='currentSpeProSeqPath_ByProteins')
    parser.add_argument('-ofa','--currentSpe_OrthologousGroup_Fa_path', type=str, help='currentSpe_OrthologousGroup_Fa_path')
    parser.add_argument('-out','--currentSpe_phmmer_outPath', type=str, help='currentSpe_phmmer_outPath')
    parser.add_argument('-log','--currentSpe_phmmer_logPath', type=str, help='currentSpe_phmmer_logPath')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    parser.add_argument('-ut','--code_utilities_folder', type=str, help='code_utilities_folder')   
    parser.add_argument('-ph','--phmmer_path', type=str, help='phmmer_path')   

    args = parser.parse_args()
    currentSpeProSeqPath_ByProteins=args.currentSpeProSeqPath_ByProteins
    currentSpe_OrthologousGroup_Fa_path=args.currentSpe_OrthologousGroup_Fa_path
    currentSpe_phmmer_outPath=args.currentSpe_phmmer_outPath
    currentSpe_phmmer_logPath=args.currentSpe_phmmer_logPath
    mp_task_nums=int(args.mp_task_nums)
    code_utilities_folder=args.code_utilities_folder
    phmmer_path=args.phmmer_path
                        
    currentSpe_OrthologousGroup_Fa_files=glob.glob(os.path.join(currentSpe_OrthologousGroup_Fa_path,"*"))
    print("len(currentSpe_OrthologousGroup_Fa_files):",len(currentSpe_OrthologousGroup_Fa_files))
    currentSpe_OrthologousGroup_Fa_idxs=[i+1 for i in range(len(currentSpe_OrthologousGroup_Fa_files))]

    
    
    currentSpe_OrthologousGroup_Fa_idxs_ArgFroPhmmer=[(code_utilities_folder,phmmer_path,idx,currentSpe_OrthologousGroup_Fa_path,currentSpeProSeqPath_ByProteins,currentSpe_phmmer_outPath,currentSpe_phmmer_logPath) for idx in currentSpe_OrthologousGroup_Fa_idxs]
    pool=mp.Pool(mp_task_nums) 
    pool.map(fun_newSingleMSA_EggNOG_OrthologousGroup_phmmer,currentSpe_OrthologousGroup_Fa_idxs_ArgFroPhmmer)
    pool.close() 


import argparse
import os 
import glob
import multiprocessing as mp


from create_singleMSA import hmmalign

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-hmm','--currentSpe_hmm_profiles_path', type=str, help='currentSpe_hmm_profiles_path')
    parser.add_argument('-a','--currentSpe_hmm_align_path', type=str, help='currentSpe_hmm_align_path')
    parser.add_argument('-ofa','--currentSpe_OrthologousGroup_Fa_path', type=str, help='currentSpe_OrthologousGroup_Fa_path')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    parser.add_argument('-align','--hmmalign_path', type=str, help='hmmalign_path')

    args = parser.parse_args()
    currentSpe_hmm_profiles_path=args.currentSpe_hmm_profiles_path
    currentSpe_hmm_align_path=args.currentSpe_hmm_align_path
    currentSpe_OrthologousGroup_Fa_path=args.currentSpe_OrthologousGroup_Fa_path
    mp_task_nums=int(args.mp_task_nums)
    hmmalign_path=args.hmmalign_path




    hmm_files=glob.glob(os.path.join(currentSpe_hmm_profiles_path,"*.hmm"))
    hmm_files=[os.path.basename(f) for f in hmm_files]
    hmm_files=[f[:-4] for f in hmm_files]

    print("len(hmm_files):",len(hmm_files),hmm_files[0])


    hmm_files_ArgForhmmalign=[(name,currentSpe_hmm_align_path,currentSpe_hmm_profiles_path,currentSpe_OrthologousGroup_Fa_path,hmmalign_path) for name in hmm_files]
    pool = mp.Pool(50)
    pool.map(hmmalign, hmm_files_ArgForhmmalign)
    pool.close()


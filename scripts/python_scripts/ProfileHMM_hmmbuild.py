
import argparse
import os 
import glob
import multiprocessing as mp


from create_singleMSA import hmmbuild

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-hmm','--currentSpe_hmm_profiles_path', type=str, help='currentSpe_hmm_profiles_path')
    parser.add_argument('-c','--currentSpe_ClustoMSA_path', type=str, help='currentSpe_ClustoMSA_path')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    parser.add_argument('-hmmp','--hmmbuild_path', type=str, help='hmmbuild_path')

    args = parser.parse_args()
    currentSpe_hmm_profiles_path=args.currentSpe_hmm_profiles_path
    currentSpe_ClustoMSA_path=args.currentSpe_ClustoMSA_path
    mp_task_nums=int(args.mp_task_nums)
    hmmbuild_path=args.hmmbuild_path
    
    
    
    # now build hmm profile use seed alignment  using next block
    ClustoMSA_files=glob.glob(os.path.join(currentSpe_ClustoMSA_path,"*"))
    ClustoMSA_files=[os.path.basename(f) for f in ClustoMSA_files]
    print("len(ClustoMSA_files):",len(ClustoMSA_files),ClustoMSA_files[0])

    ClustoMSA_files_ArgForhmmbuild=[(name,currentSpe_hmm_profiles_path,currentSpe_hmm_profiles_path,currentSpe_ClustoMSA_path,hmmbuild_path) for name in ClustoMSA_files]
    pool = mp.Pool(50)
    pool.map(hmmbuild, ClustoMSA_files_ArgForhmmbuild)
    pool.close()


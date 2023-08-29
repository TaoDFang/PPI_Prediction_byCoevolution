
import argparse
import os 
import glob
import multiprocessing as mp


from create_singleMSA import ClustoMSA

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-pf','--currentSpe_phmmer_OrthologousGroup_path', type=str, help='currentSpe_phmmer_OrthologousGroup_path')
    parser.add_argument('-c','--currentSpe_ClustoMSA_path', type=str, help='currentSpe_ClustoMSA_path')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    parser.add_argument('-cp','--clustalo_path', type=str, help='clustalo_path')

    args = parser.parse_args()
    currentSpe_phmmer_OrthologousGroup_path=args.currentSpe_phmmer_OrthologousGroup_path
    currentSpe_ClustoMSA_path=args.currentSpe_ClustoMSA_path
    mp_task_nums=int(args.mp_task_nums)
    clustalo_path=args.clustalo_path

       
    phmmer_OrthologousGroup_names=glob.glob(currentSpe_phmmer_OrthologousGroup_path+"*.fa")
    phmmer_OrthologousGroup_names=[os.path.basename(f) for f in phmmer_OrthologousGroup_names]
    phmmer_OrthologousGroup_names=[f[:-3] for f in phmmer_OrthologousGroup_names]

    print("len(phmmer_OrthologousGroup_names),phmmer_OrthologousGroup_names[0]:",len(phmmer_OrthologousGroup_names),phmmer_OrthologousGroup_names[0])


    phmmer_OrthologousGroup_names_ArgFroClustoMSA=[(name,currentSpe_phmmer_OrthologousGroup_path,currentSpe_ClustoMSA_path,clustalo_path) for name in phmmer_OrthologousGroup_names]
    pool = mp.Pool(mp_task_nums)
    pool.map(ClustoMSA, phmmer_OrthologousGroup_names_ArgFroClustoMSA)
    pool.close()



# CPU times: user 34.3 s, sys: 1min 34s, total: 2min 9s
# Wall time: 1h 11min 7s

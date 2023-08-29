import argparse
import os 
import glob
import pandas as pd 
import multiprocessing as mp
from multiprocessing import get_context


from create_pairedMSA import  getSameProteinRatio


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-nf90csv','--pairedMSA_Nf90_csv', type=str, help='pairedMSA_Nf90_csv')
    parser.add_argument('-nf90f','--pairedMSA_Nf90_folder', type=str, help='pairedMSA_Nf90_folder')
    parser.add_argument('-scsv','--pairedMSA_sameProteinRatio_csv', type=str, help='pairedMSA_sameProteinRatio_csv')
    parser.add_argument('-acsv','--PPIInfoBeforeCoEvoComp_csv', type=str, help='PPIInfoBeforeCoEvoComp_csv')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')
    
    
    args = parser.parse_args()
    pairedMSA_Nf90_csv=args.pairedMSA_Nf90_csv
    pairedMSA_Nf90_folder=args.pairedMSA_Nf90_folder
    pairedMSA_sameProteinRatio_csv=args.pairedMSA_sameProteinRatio_csv
    PPIInfoBeforeCoEvoComp_csv=args.PPIInfoBeforeCoEvoComp_csv
    mp_task_nums=int(args.mp_task_nums)
    
    # get postivte and negative ppi  second setp (after get all needed final paired MSA)
    # by  same protein ratio on two side of MSA 


    # first need to update postive and negative ppi to only use pp that have large Nf90 value 
    pairedMSA_Nf90_frame = pd.read_csv(pairedMSA_Nf90_csv,header=None,index_col=None,sep="\t")
    print("pairedMSA_Nf90_frame.shape:",pairedMSA_Nf90_frame.shape)

    pairedMSA_Nf90_list=pairedMSA_Nf90_frame.iloc[:,[0,1]].values.tolist()
    currentSpe_allPPIs_pps=[(p1,p2) for p1, p2 in pairedMSA_Nf90_list]
    print("len(currentSpe_allPPIs_pps):",len(currentSpe_allPPIs_pps))



    #heck paired MSA of homologus , if same proteins from same species in both side
    if not os.path.exists(pairedMSA_sameProteinRatio_csv):
        with open(pairedMSA_sameProteinRatio_csv, 'w') as fp:
            pass  
    #in the nextflow piple, this block is meaningless as everytime the process start from begining and preprall all possible ppi 
    pairedMSA_sameProteinRatio_dict=dict()
    # if os.path.getsize(pairedMSA_sameProteinRatio_csv) == 0:
    #     pairedMSA_sameProteinRatio_dict=dict()
    # else:
    #     pairedMSA_sameProteinRatio_frame = pd.read_csv(pairedMSA_sameProteinRatio_csv,header=None,index_col=None,sep="\t")
    #     pairedMSA_sameProteinRatio_list=pairedMSA_sameProteinRatio_frame.values.tolist()
    #     pairedMSA_sameProteinRatio_dict=dict([((p1,p2),r) for p1 , p2 , r in pairedMSA_sameProteinRatio_list])



    currentSpe_allPPIs_ArgForSamePro=[(p1,p2,pairedMSA_Nf90_folder) for p1, p2, in currentSpe_allPPIs_pps if (p1, p2) not in pairedMSA_sameProteinRatio_dict]
    print("len(currentSpe_allPPIs_ArgForSamePro):",len(currentSpe_allPPIs_ArgForSamePro))

    if len(currentSpe_allPPIs_ArgForSamePro)>0:
        pool=mp.Pool(mp_task_nums)
        sameProtein_ratio_results=pool.map(getSameProteinRatio,currentSpe_allPPIs_ArgForSamePro)
        pool.close() 
        sameProtein_ratio_frame=pd.DataFrame(sameProtein_ratio_results,columns=["protein1","protein2","saemProtein_ratio"])
        sameProtein_ratio_frame.to_csv(pairedMSA_sameProteinRatio_csv,mode="a",
                            header=None,index=None,sep="\t")
        print("sameProtein_ratio_frame.shape:",sameProtein_ratio_frame.shape)



    # remove protein pair with same prortein on two side of MSA 

    sameProtein_ratio_frame=pd.read_csv(pairedMSA_sameProteinRatio_csv,
                            header=None,index_col=None,sep="\t")
    print("sameProtein_ratio_frame.shape",sameProtein_ratio_frame.shape)
    large_sameProtein_ratio_frame=sameProtein_ratio_frame.loc[sameProtein_ratio_frame.iloc[:,2]>0,:]
    print("large_sameProtein_ratio_frame.shape",large_sameProtein_ratio_frame.shape)

    large_sameProtein_ratio_dict=dict([((p1,p2),r) for p1,p2,r in large_sameProtein_ratio_frame.values.tolist()])


    currentSpe_allPPIs_pps=[pp for pp in currentSpe_allPPIs_pps if pp not in large_sameProtein_ratio_dict]
    print("len(currentSpe_allPPIs_pps):",len(currentSpe_allPPIs_pps))



    currentSpe_allPPIs_beforeCoEvoComp_info=[(p1,p2, "P") for p1, p2 in currentSpe_allPPIs_pps]

    currentSpe_allPPIs_beforeCoEvoComp_frame=pd.DataFrame(currentSpe_allPPIs_beforeCoEvoComp_info,
                                                                               columns=["STRING_ID1","STRING_ID2","benchmark_status"])
    currentSpe_allPPIs_beforeCoEvoComp_frame.to_csv(PPIInfoBeforeCoEvoComp_csv,
                                                                         header=True, index=False,sep="\t")
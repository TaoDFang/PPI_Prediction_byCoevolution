import argparse
import os 
import glob
import multiprocessing as mp


from IntergrateBestHomologousPPCoEvo_unNameSorted import getEggNOG_homologousNameUnSortedPP_COG2PP_reversedPP_dict


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('-s','--spe_list', type=str, help='spe_list')
    parser.add_argument('-egg','--EggNOG_groupPath', type=str, help='EggNOG_groupPath')
    parser.add_argument('-t','--homologous_COG2PP_path', type=str, help='homologous_COG2PP_path')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')

    args = parser.parse_args()
    spe_list=args.spe_list
    EggNOG_groupPath=args.EggNOG_groupPath
    homologous_COG2PP_path=args.homologous_COG2PP_path
    mp_task_nums=int(args.mp_task_nums)



    max_level="2"
    Spe_For_getEggNOG_homologousNameUnSortedPP_COG2PP_reversedPP_dict=spe_list #['511145','1274374','105422','411476']
    Arg_For_getEggNOG_homologousNameUnSortedPP_COG2PP_reversedPP_dict=[(spe,max_level,EggNOG_groupPath,homologous_COG2PP_path) for spe in Spe_For_getEggNOG_homologousNameUnSortedPP_COG2PP_reversedPP_dict]
    
    mp_task_nums=max(mp_task_nums,len(spe_list))
    pool = mp.Pool(mp_task_nums)
    pool.map(getEggNOG_homologousNameUnSortedPP_COG2PP_reversedPP_dict, Arg_For_getEggNOG_homologousNameUnSortedPP_COG2PP_reversedPP_dict)
    pool.close()


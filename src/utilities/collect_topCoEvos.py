import os 
import numpy as np
import pandas as pd
from Bio import AlignIO
import multiprocessing as mp 
import sys
import shutil
from pathlib import Path

def get_pydcaFNAPC_betArray(Coevo_path,Query_pro1,Query_pro2,
                            L1,suffix):

    data_fileName=Coevo_path+Query_pro1+"and"+Query_pro2+suffix+".npz"
    #data_fileName=data_array_folder+Query_pro1+"and"+Query_pro2+".npy"
    #print(data_fileName)
    try:
        data_array_Dig=np.load(data_fileName)['arr_0']
    #except Exception as e
    except ValueError as e :   # deal with multiple error 
        
#         #this is to for debugging , comment them for collectin data
#         #ValueError('Cannot load file containing pickled data when allow_pickle=False'), get this by return e ,show it 
#         parent_folder=Path(Coevo_path).parent.absolute()
#         corrupted_DCA_coevolutoin_path=os.path.join(parent_folder,"corrupted_coevolutoin_result_DCA/")
#         corrupted_data_fileName=corrupted_DCA_coevolutoin_path+Query_pro1+"and"+Query_pro2+suffix+".npz"
#         dest = shutil.move(data_fileName, corrupted_data_fileName)      
# #         print(f"excepttion error:{e}")
# #         print(f"corrupted_data_fileName:{corrupted_data_fileName}")
# #         print(f"dest:{dest}")
        print(data_fileName)
        return None # for now just to check how many corrupt data i have 
    except IOError as e:
        print(data_fileName)
        return None


    data_array =data_array_Dig.T + data_array_Dig
    np.fill_diagonal(data_array,0)
    bet_data_array=data_array[:L1,L1:]
    
    return bet_data_array




def get_topRankingBetValue_dict(record,topNum=50):
    Query_pro1, Query_pro2,L1,L2,CoEvo_path,suffix=record
    #print(Query_pro1, Query_pro2,L1,L2)
    
    if suffix =="_pydcaFNAPC_array":
        bet_data_array=get_pydcaFNAPC_betArray(CoEvo_path,Query_pro1,Query_pro2,
                            L1,suffix)
    if suffix =="_apc_allResidues":
        bet_data_array=get_pydcaFNAPC_betArray(CoEvo_path,Query_pro1,Query_pro2,
                            L1,suffix)
        
    if suffix=="_alphafold_prob12":
        bet_data_array=get_alphafoldPro12_betArray(CoEvo_path,Query_pro1,Query_pro2,
                            L1)
    if suffix=="_alphafold_inverseMinDist":
        bet_data_array=get_alphafoldinverseMinDist_betArray(CoEvo_path,Query_pro1,Query_pro2,
                            L1)
    if bet_data_array is None: 
        return None

    ascending_bet_data_array=np.sort(bet_data_array.flatten())
    descending_bet_data_array=ascending_bet_data_array[::-1]


    returnList=[Query_pro1,Query_pro2]
    #here there one problem. if there are same top_value
    # for each top value, following for loop will find all position of this top value, cause many repeats
    # so here use set to make use each top vop_value only used for once 
    esisted_top_values=set()
    for j in range(topNum):
        top_value=descending_bet_data_array[j]
        if top_value not in esisted_top_values:
            topIdx=np.where(bet_data_array==top_value)
            for i in range(len(topIdx[0])):
                topRow,topCol=topIdx[0][i],topIdx[1][i]
                topCol += L1
                #print(top_value,topRow,topCol)
                if len(returnList)<(2+topNum*3):
                    returnList.extend((top_value,topRow,topCol))
            esisted_top_values.add(top_value)

        if len(returnList)>(2+topNum*3):
            break
            
    return(returnList)


def get_topRanking_CoEvo_file(topRanking_CoEvo_file,coevolutoin_path,coevo_suffix,allPPI_info, 
                              returnDic=False,overwrite=True,use_multiprocessing=True):
    # here use_multiprocessing coulde be bool(True or False ), or integeger number which represent how many cpus planed to use 
    if  overwrite:
        print("topRanking_CoEvo_file:",topRanking_CoEvo_file)


        allPPI_info_ArgForTops=[[p1,p2,L1,L2,coevolutoin_path,coevo_suffix] for p1,p2,L1,L2 in allPPI_info]
    
        if use_multiprocessing:
            pool=mp.Pool(use_multiprocessing)
            top_CoEvo_list=pool.map(get_topRankingBetValue_dict,allPPI_info_ArgForTops)
            pool.close() 
        else:
            top_CoEvo_list=[get_topRankingBetValue_dict(l) for l in allPPI_info_ArgForTops]
        
        print(f"len(top_CoEvo_list):{len(top_CoEvo_list)}")
        top_CoEvo_list=[l for l in top_CoEvo_list if l is not None]
        
        print(f"len(top_CoEvo_list):{len(top_CoEvo_list)}")
        print(sys.getsizeof(top_CoEvo_list)/(1024*1024*1024))
        
        top_CoEvo_frame=pd.DataFrame(top_CoEvo_list)
        top_CoEvo_frame.to_csv(topRanking_CoEvo_file,header=None, index=None,sep="\t")


    top_CoEvo_frame=pd.read_csv(topRanking_CoEvo_file,header=None, index_col=None,sep="\t")

    print("top_CoEvo_frame.shape:",top_CoEvo_frame.shape)
    
    if returnDic:
        top_CoEvo_list=top_CoEvo_frame.values.tolist()
        top_CoEvo_dict={(l[0],l[1]): l[2:] for l in top_CoEvo_list} 
        print(f"sys.getsizeof(top_CoEvo_list)/(1024*1024*1024):{sys.getsizeof(top_CoEvo_list)/(1024*1024*1024)}")
        return (top_CoEvo_dict)
    else:
        print(f"sys.getsizeof(top_CoEvo_frame)/(1024*1024*1024):{sys.getsizeof(top_CoEvo_frame)/(1024*1024*1024)}")
        return(top_CoEvo_frame)
    
    
    
    
    
def collect_topCoEvos_OnlyTopPosNeg(topCoEvoScore_frame,allPPI_allInfoFrame,CoEvo_type,DCA_thres,
                                    CoEvo_number,selCoEvo_number,ML_inputPath,
                                    ifReCollect=False):
    
    X_fileName=ML_inputPath+"DCA_thres_"+str(DCA_thres)+"_NonPara_X_top"+CoEvo_type+"s_Num"+str(CoEvo_number)+".npy"
    Y_fileName=ML_inputPath+"DCA_thres_"+str(DCA_thres)+"_NonPara_Y_top"+CoEvo_type+"s_Num"+str(CoEvo_number)+".npy"
    if ifReCollect:
        print("Re collect data ")
        print("X_fileName:",X_fileName)
        print("Y_fileName:",Y_fileName)
        CoEvo_index=[3*n for n in range(CoEvo_number)]

        top_CoEvo_list=topCoEvoScore_frame.values.tolist()
        top_CoEvo_list_tuple=[((l[0:2]),(l[2:]))for l in top_CoEvo_list]
        top_CoEvo_list_tuple=[(tuple(l1),tuple(l2)) for l1,l2 in top_CoEvo_list_tuple]
        top_CoEvo_dict=dict(top_CoEvo_list_tuple)


        NonPara_allPPI_info=allPPI_allInfoFrame.loc[:,["STRING_ID1","STRING_ID2","benchmark_status"]].values.tolist()



        OnlyTopPosNeg_NonPara_X=np.zeros((len(NonPara_allPPI_info),CoEvo_number))
        OnlyTopPosNeg_NonPara_Y=np.zeros((len(NonPara_allPPI_info),1))

        for idx,l in enumerate(NonPara_allPPI_info):
            p1, p2 , label=l
            numeric_lable=1 if label=="P" else 0

            cur_coevos=top_CoEvo_dict[(p1,p2)]
            OnlyTopPosNeg_NonPara_X[idx,:]=[cur_coevos[i] for i in CoEvo_index]
            OnlyTopPosNeg_NonPara_Y[idx]=numeric_lable


        np.save(X_fileName,OnlyTopPosNeg_NonPara_X)
        np.save(Y_fileName,OnlyTopPosNeg_NonPara_Y)

    OnlyTopPosNeg_NonPara_XtopCoEvos=np.load(X_fileName)
    OnlyTopPosNeg_NonPara_XtopCoEvos=OnlyTopPosNeg_NonPara_XtopCoEvos[:,0:selCoEvo_number]

    OnlyTopPosNeg_NonPara_YtopCoEvos=np.load(Y_fileName)
    OnlyTopPosNeg_NonPara_YtopCoEvos=np.reshape(OnlyTopPosNeg_NonPara_YtopCoEvos,(OnlyTopPosNeg_NonPara_YtopCoEvos.shape[0],))#

    return(OnlyTopPosNeg_NonPara_XtopCoEvos,OnlyTopPosNeg_NonPara_YtopCoEvos)





def get_topRankingBetValue_dict_PosInSingleMSA_framArray(bet_data_array,topNum=50):

    ascending_bet_data_array=np.sort(bet_data_array.flatten())
    descending_bet_data_array=ascending_bet_data_array[::-1]
    returnList=[]
    esisted_top_values=set()
    for j in range(topNum):
        top_value=descending_bet_data_array[j]
        if top_value not in esisted_top_values:
            topIdx=np.where(bet_data_array==top_value)
            for i in range(len(topIdx[0])):
                topRow,topCol=topIdx[0][i],topIdx[1][i]
                if len(returnList)<(2+topNum*3):
                    returnList.extend((top_value,topRow,topCol))
            esisted_top_values.add(top_value)

        if len(returnList)>(2+topNum*3):
            break
            
    return(returnList)
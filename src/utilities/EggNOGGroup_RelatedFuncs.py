import sys
import os 
import pandas as pd
import multiprocessing as mp 
import numpy as np 

from collections import defaultdict
import copy

from sklearn.model_selection import train_test_split


def flat2list(ll1,ll2):
    l1=copy.deepcopy(ll1)
    l2=copy.deepcopy(ll2)
    
    if (len(l1)==0) or (len(l2)==0):
        print("something wrong betweenn protein to cog group mapping ")
    l1.extend(l2) # here actully modificion list outside of function when they are input arguments 
    # https://stackoverflow.com/questions/31435603/python-modify-global-list-inside-a-function
    return(l1)

def return_cogPairs_fromproPair(EggNOG_group_level2,currentSpe_TaxID,allPPI_allInfo_frame):
    EggNOG_group_level2_speID=EggNOG_group_level2.iloc[:,2].values.tolist()
    EggNOG_group_level2_speID=[str(sid) for sid in EggNOG_group_level2_speID]
    
    CurrentSpe_EggNOG_group_level2_idx=[idx for idx,sid in  enumerate(EggNOG_group_level2_speID) if (sid==currentSpe_TaxID)]
    CurrentSpe_EggNOG_group_level2=EggNOG_group_level2.iloc[CurrentSpe_EggNOG_group_level2_idx,:]
    #CurrentSpe_EggNOG_group_level2=CurrentSpe_EggNOG_group_level2.sort_values(by=1)



    CurrentSpe_EggNOG_group_level2_values=CurrentSpe_EggNOG_group_level2.iloc[:,[1,3]].values.tolist()
    CurrentSpe_EggNOG_group_level2_pro2cog=defaultdict(list)

    for cog, pro  in CurrentSpe_EggNOG_group_level2_values:
        CurrentSpe_EggNOG_group_level2_pro2cog[pro].append(cog)

    print("len(CurrentSpe_EggNOG_group_level2_pro2cog):",len(CurrentSpe_EggNOG_group_level2_pro2cog))


    allPPI_pps=allPPI_allInfo_frame.loc[:,["STRING_ID1","STRING_ID2"]].values.tolist()
    allPPI_cogs=[flat2list(CurrentSpe_EggNOG_group_level2_pro2cog[p1][:],CurrentSpe_EggNOG_group_level2_pro2cog[p2][:]) for p1, p2 in allPPI_pps]

    print("len(allPPI_cogs):",len(allPPI_cogs) )
    return(allPPI_cogs)


def sepCogPairs_train_test_split_getTrainandTestPPTuples(allPPI_cogs,all_PPTuples,test_size=0.20, random_state=0):
    cogs_train,cogs_test,all_PPTuples_train,all_PPTuples_test = train_test_split(allPPI_cogs,all_PPTuples, test_size=test_size, random_state=random_state)

    
    print("type(cogs_train),type(cogs_test):",type(cogs_train),type(cogs_test))
    
    print(cogs_train[0:3])
    cogs_train_dict={tuple(cog):1 for cog in cogs_train}
    print("len(cogs_train_dict):",len(cogs_train_dict))
    Moved_test_idx=[]
    Kepted_test_idx=[]
    for i,cogs in enumerate(cogs_test):
        if tuple(cogs) in cogs_train_dict:
            Moved_test_idx.append(i)
        else:
            Kepted_test_idx.append(i)
    print("len(cogs_test),len(Moved_test_idx),len(Kepted_test_idx):", len(cogs_test),len(Moved_test_idx),len(Kepted_test_idx))

    
    all_PPTuples_train=all_PPTuples_train+[all_PPTuples_test[i] for i in Moved_test_idx]
    all_PPTuples_test=[all_PPTuples_test[i] for i in Kepted_test_idx]
    
    print("len(all_PPTuples_train),len(all_PPTuples_test):",len(all_PPTuples_train),len(all_PPTuples_test))
    return(all_PPTuples_train,all_PPTuples_test)

def sepCogPairs_train_test_split(allPPI_cogs,OnlyTopPosNeg_NonPara_XtopDCAs, OnlyTopPosNeg_NonPara_YtopDCAs,test_size=0.20, random_state=0):
    
    print("OnlyTopPosNeg_NonPara_XtopDCAs.ndim:",OnlyTopPosNeg_NonPara_XtopDCAs.ndim)
    if OnlyTopPosNeg_NonPara_XtopDCAs.ndim==1:
        OnlyTopPosNeg_NonPara_XtopDCAs=np.reshape(OnlyTopPosNeg_NonPara_XtopDCAs,(-1,1))
        
    cogs_train,cogs_test,XtopDCAs_train, XtopDCAs_test, ytopDCAs_train, ytopDCAs_test = train_test_split(allPPI_cogs,OnlyTopPosNeg_NonPara_XtopDCAs, OnlyTopPosNeg_NonPara_YtopDCAs, test_size=test_size, random_state=random_state)
    print("XtopDCAs_train.shape,XtopDCAs_test.shape,ytopDCAs_train.shape,ytopDCAs_test.shape:",XtopDCAs_train.shape,XtopDCAs_test.shape,ytopDCAs_train.shape,ytopDCAs_test.shape)

    print(cogs_train[0:3])
    cogs_train_dict={tuple(cog):1 for cog in cogs_train}
    print("len(cogs_train_dict):",len(cogs_train_dict))
    Moved_test_idx=[]
    Kepted_test_idx=[]
    for i,cogs in enumerate(cogs_test):
        if tuple(cogs) in cogs_train_dict:
            Moved_test_idx.append(i)
        else:
            Kepted_test_idx.append(i)
    print("len(cogs_test),len(Moved_test_idx),len(Kepted_test_idx):", len(cogs_test),len(Moved_test_idx),len(Kepted_test_idx))

    
    XtopDCAs_train=np.concatenate([XtopDCAs_train,XtopDCAs_test[Moved_test_idx,:]],axis=0)
    ytopDCAs_train=np.concatenate([ytopDCAs_train,ytopDCAs_test[Moved_test_idx]],axis=0)
    XtopDCAs_test=XtopDCAs_test[Kepted_test_idx,:]
    ytopDCAs_test=ytopDCAs_test[Kepted_test_idx]
    
    cogs_train_array=np.array(cogs_train,dtype=object)
    cogs_test_array=np.array(cogs_test,dtype=object)
    cogs_train_array=np.concatenate([cogs_train_array,cogs_test_array[Moved_test_idx]],axis=0)
    
    
    return(cogs_train_array,XtopDCAs_train, XtopDCAs_test, ytopDCAs_train, ytopDCAs_test)
    
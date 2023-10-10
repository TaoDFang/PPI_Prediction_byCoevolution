import pandas as pd 
import numpy as np
import copy 
import random
from collections import defaultdict



def removeTrainPPfrom_allPredictionResults_ByRemovedPosRatio(result_list,pp_list,trainPP_dict,allPos_dict,thres,):
    def returnListByIdx(input_list,idx_list):
        return([input_list[idx] for idx in idx_list])
    # here notice pp names  in pp_list and trainPP_dict should be both sorted 
    assert len(result_list)==len(pp_list)
    print(len(result_list),len(pp_list))
    
    # sep results_list to pos and neg using benchmakr allPosDict
    result_posIdx=[idx for idx,pp in enumerate(pp_list) if pp in allPos_dict ]
    result_negIdx=[idx for idx,pp in enumerate(pp_list) if pp not in allPos_dict]
    print("len(result_posIdx),len(result_negIdx):",len(result_posIdx),len(result_negIdx))
    
    result_list_pos=returnListByIdx(result_list,result_posIdx)
    pp_list_pos=returnListByIdx(pp_list,result_posIdx)
    
    result_list_neg=returnListByIdx(result_list,result_negIdx)
    pp_list_neg=returnListByIdx(pp_list,result_negIdx)
    

    
    # now remove training propertionaly 
    kept_posIdx=[idx for idx,pp in enumerate(pp_list_pos) if pp not in trainPP_dict]
    result_list_pos_kept=returnListByIdx(result_list_pos,kept_posIdx)
    pp_list_pos_kept=returnListByIdx(pp_list_pos,kept_posIdx)
    print("len(pp_list_pos_kept):",len(pp_list_pos_kept))
    
    kept_negIdx=[idx for idx,pp in enumerate(pp_list_neg) if pp not in trainPP_dict]
    result_list_neg_kept=returnListByIdx(result_list_neg,kept_negIdx)
    pp_list_neg_kept=returnListByIdx(pp_list_neg,kept_negIdx)
    print("len(pp_list_neg_kept):,before fix ratio",len(pp_list_neg_kept))
    
    nonOverlaped_pos_len=len(kept_posIdx)
    kept_ratio=nonOverlaped_pos_len/len(pp_list_pos)
    print("kept_ratio:",kept_ratio)
    
    kepted_negIdx_num=int(len(result_list_neg)*(kept_ratio))#-len(removed_negIdx)
    
    shuffle_negIdx=list(range(len(pp_list_neg_kept)))
    random.Random(10).shuffle(shuffle_negIdx)
    kept_negIdx=shuffle_negIdx[0:kepted_negIdx_num]
    result_list_neg_kept=returnListByIdx(result_list_neg_kept,kept_negIdx)
    pp_list_neg_kept=returnListByIdx(pp_list_neg_kept,kept_negIdx)
    
    print("len(result_list_pos_kept),len(result_list_neg_kept):",len(result_list_pos_kept),len(result_list_neg_kept))
    
    result_list_kept=result_list_pos_kept+result_list_neg_kept
    pp_list_kept=pp_list_pos_kept+pp_list_neg_kept
    print(len(result_list_kept),len(pp_list_kept))
    
    if thres is not None: 
        print("futher filtering by some threshold")

        final_idx=[idx for idx,v in enumerate(result_list_kept) if v>thres]
        final_result_pos=[result_list_kept[idx] for idx in final_idx]  # ie, final prediction                  
        final_pp_pos=[pp_list_kept[idx] for idx in final_idx]                

        print(len(final_result_pos),len(final_pp_pos))
        return(final_result_pos,final_pp_pos)
    else:
        return(result_list_kept,pp_list_kept)
    
    
    
def phylaIntegration_replacingDCAScores(ori_inputMatrix,fillDCAValue,fillMissingValue,):

    inputMatrix=copy.deepcopy(ori_inputMatrix)
    inputMatrix[np.where(inputMatrix!=fillMissingValue)] = fillDCAValue
    return(inputMatrix)

def phylaIntegration_replacingSubjectSpeDCAScores(ori_inputMatrix,fillDCAValue,fillMissingValue,SubjectFeaNum=1):

    inputMatrix=copy.deepcopy(ori_inputMatrix)
    inputMatrix[:,0:SubjectFeaNum]= fillDCAValue
    return(inputMatrix)

def phylaIntegration_replacingOtherPhalaDCAScores(ori_inputMatrix,fillDCAValue,fillMissingValue,SubjectFeaNum=1):

    inputMatrix=copy.deepcopy(ori_inputMatrix)
    inputMatrix_col0 = copy.deepcopy(inputMatrix[:,0:SubjectFeaNum])
    
    inputMatrix[np.where(inputMatrix!=fillMissingValue)] = fillDCAValue
    inputMatrix[:,0:SubjectFeaNum]= inputMatrix_col0
    return(inputMatrix)



def getAllIdxOf_fixed_negVSpos_ratio_keepCombination(ori_XMatrix,ori_YMatrix,fixedCombinatioNegVSposRratio):
    XMatrix=copy.deepcopy(ori_XMatrix)
    YMatrix=copy.deepcopy(ori_YMatrix)

    posIdx_dict=defaultdict(list)
    negIdx_dict=defaultdict(list)

    for i in range(0,XMatrix.shape[0]):
        i_phyla=XMatrix[i,:]
        i_phyla_label=[str(p) for p in i_phyla]
        i_phyla_label="_".join(i_phyla_label)

        i_y=YMatrix[i]
        if i_y==1:
            posIdx_dict[i_phyla_label].append(i)
        elif i_y==0:
            negIdx_dict[i_phyla_label].append(i)
    
    print("posIdx_dict")
    for k, v in posIdx_dict.items():
        print(k,len(v),v[0:10])
        
    print("negIdx_dict")
    for k, v in negIdx_dict.items():
        print(k,len(v),v[0:10])
        
    for k, v in negIdx_dict.items():
        random.Random(10).shuffle(v)
            
    print("After shuffling negIdx_dict")
    for k, v in negIdx_dict.items():
        print(k,len(v),v[0:10])
        
    revised_negIdx_dict=dict()
    for k, v in negIdx_dict.items():
        revised_negNum=int(fixedCombinatioNegVSposRratio*len(posIdx_dict[k]))
        print("revised_negNum:",k,len(posIdx_dict[k]),revised_negNum)
        revised_negIdx_dict[k]=v[0:revised_negNum]
    
    print("After fixedNegVSposRratio negIdx_dict")
    for k, v in revised_negIdx_dict.items():
        print(k,len(v),v[0:10])
        
    final_returnIdx=list()
    for k, v in posIdx_dict.items():
        final_returnIdx.extend(v)
        
    for k, v in revised_negIdx_dict.items():
        final_returnIdx.extend(v)
        
    final_returnIdx=sorted(final_returnIdx)
    
    print("len(final_returnIdx)):",len(final_returnIdx))
    return(final_returnIdx)


def getAllIdxOf_fixed_negVSpos_ratio_keepCombination_from_OriginalFrame(XBestHomologousDCAs_original,
                                                                        YBestHomologousDCAs_original,
                                                                        fillMissingValue=-1,
                                                                       ifplot=False):
    XBestHomologousDCAs_replace=copy.deepcopy(XBestHomologousDCAs_original)
    YBestHomologousDCAs_replace=copy.deepcopy(YBestHomologousDCAs_original)
    print("YBestHomologousDCAs_replace.shape:",YBestHomologousDCAs_replace.shape)

    XBestHomologousDCAs_replace[XBestHomologousDCAs_replace!=fillMissingValue]=1

    beforeFixRatio_test_Phaly_posPPInum=defaultdict(int)
    beforeFixRatio_test_Phaly_negPPInum=defaultdict(int)

    beforeFixRatio_test_PhalyCombination_posPPInum =defaultdict(int)
    beforeFixRatio_test_PhalyCombination_negPPInum=defaultdict(int)

    for i in range(0,XBestHomologousDCAs_replace.shape[0]):
        i_phyla=XBestHomologousDCAs_replace[i,:]
        i_phyla_ratio=sum(i_phyla==fillMissingValue)/XBestHomologousDCAs_replace.shape[1]
        i_phyla_label=[str(p) for p in i_phyla]
        i_phyla_label="_".join(i_phyla_label)

        i_y=YBestHomologousDCAs_replace[i]
        if i_y==1:
            beforeFixRatio_test_Phaly_posPPInum[i_phyla_ratio] +=1
            beforeFixRatio_test_PhalyCombination_posPPInum[i_phyla_label] +=1
        elif i_y==0:
            beforeFixRatio_test_Phaly_negPPInum[i_phyla_ratio] +=1
            beforeFixRatio_test_PhalyCombination_negPPInum[i_phyla_label] +=1


    beforeFixRatio_test_Phaly_negVSpos_ratio={minus1_Ratio:number/(beforeFixRatio_test_Phaly_posPPInum[minus1_Ratio]) for minus1_Ratio,number in beforeFixRatio_test_Phaly_negPPInum.items()}
    print("beforeFixRatio_test_Phaly_posPPInum, beforeFixRatio_test_Phaly_negPPInum,beforeFixRatio_test_Phaly_negVSpos_ratio:",beforeFixRatio_test_Phaly_posPPInum, beforeFixRatio_test_Phaly_negPPInum,beforeFixRatio_test_Phaly_negVSpos_ratio)



    beforeFixRatio_test_Phaly_PPInum_list=list()
    for key in sorted(list(beforeFixRatio_test_Phaly_posPPInum.keys())):
        print(key)
        beforeFixRatio_test_Phaly_PPInum_list.append([key,beforeFixRatio_test_Phaly_posPPInum[key],"P"])
        beforeFixRatio_test_Phaly_PPInum_list.append([key,beforeFixRatio_test_Phaly_negPPInum[key],"N"])
    beforeFixRatio_test_Phaly_PPInum_frame=pd.DataFrame(beforeFixRatio_test_Phaly_PPInum_list,columns=["-1Ratio","number","benchmark_status"])



    #****************revise ratio at each -1 group , and consider order and combination of -1******************
    beforeFixRatio_test_PhalyCombination_negVSpos_ratio={minus1_Ratio:number/(beforeFixRatio_test_PhalyCombination_posPPInum[minus1_Ratio]) for minus1_Ratio,number in beforeFixRatio_test_PhalyCombination_negPPInum.items()}
    for  key in beforeFixRatio_test_PhalyCombination_posPPInum.keys():
        print(key,beforeFixRatio_test_PhalyCombination_posPPInum[key], beforeFixRatio_test_PhalyCombination_negPPInum[key],beforeFixRatio_test_PhalyCombination_negVSpos_ratio[key])


    min_fixedCombination_negVSpos_ratio=beforeFixRatio_test_PhalyCombination_negVSpos_ratio[min(beforeFixRatio_test_PhalyCombination_negVSpos_ratio,key=beforeFixRatio_test_PhalyCombination_negVSpos_ratio.get)]
    print("min_fixedCombination_negVSpos_ratio:",min_fixedCombination_negVSpos_ratio)



    fixedNegVSposRratio_keepCombination_allIdx=getAllIdxOf_fixed_negVSpos_ratio_keepCombination(ori_XMatrix=XBestHomologousDCAs_replace,
                                     ori_YMatrix=YBestHomologousDCAs_replace,
                                     fixedCombinatioNegVSposRratio=min_fixedCombination_negVSpos_ratio,
                                                                                               )

    # here show for visulisation reason 
    XMarix_fixedNegVSposRratio_keepCombination=XBestHomologousDCAs_replace[fixedNegVSposRratio_keepCombination_allIdx,:]
    YMarix_fixedNegVSposRratio_keepCombination=YBestHomologousDCAs_replace[fixedNegVSposRratio_keepCombination_allIdx]
    print("YMarix_fixedNegVSposRratio_keepCombination.shape:",YMarix_fixedNegVSposRratio_keepCombination.shape)

    test_Phaly_posPPInum_fixedNegVSposRratio_keepCombination =defaultdict(int)
    test_Phaly_negPPInum_fixedNegVSposRratio_keepCombination=defaultdict(int)

    test_PhalyCombination_posPPInum_fixedNegVSposRratio_keepCombination =defaultdict(int)
    test_PhalyCombination_negPPInum_fixedNegVSposRratio_keepCombination=defaultdict(int)

    for i in range(0,XMarix_fixedNegVSposRratio_keepCombination.shape[0]):
        i_phyla=XMarix_fixedNegVSposRratio_keepCombination[i,:]
        i_phyla_ratio=sum(i_phyla==-1)/XMarix_fixedNegVSposRratio_keepCombination.shape[1]
        i_phyla_label=[str(p) for p in i_phyla]
        i_phyla_label="_".join(i_phyla_label)

        i_y=YMarix_fixedNegVSposRratio_keepCombination[i]
        if i_y==1:
            test_Phaly_posPPInum_fixedNegVSposRratio_keepCombination[i_phyla_ratio] +=1
            test_PhalyCombination_posPPInum_fixedNegVSposRratio_keepCombination[i_phyla_label] +=1
        elif i_y==0:
            test_Phaly_negPPInum_fixedNegVSposRratio_keepCombination[i_phyla_ratio] +=1
            test_PhalyCombination_negPPInum_fixedNegVSposRratio_keepCombination[i_phyla_label] +=1

    test_PhalyCombination_negVSpos_ratio_fixedNegVSposRratio_keepCombination={minus1_Ratio:number/(test_PhalyCombination_posPPInum_fixedNegVSposRratio_keepCombination[minus1_Ratio]) for minus1_Ratio,number in test_PhalyCombination_negPPInum_fixedNegVSposRratio_keepCombination.items()}
    for  key in test_PhalyCombination_posPPInum_fixedNegVSposRratio_keepCombination.keys():
        print(key,test_PhalyCombination_posPPInum_fixedNegVSposRratio_keepCombination[key], test_PhalyCombination_negPPInum_fixedNegVSposRratio_keepCombination[key],test_PhalyCombination_negVSpos_ratio_fixedNegVSposRratio_keepCombination[key])


    test_Phaly_PPInum_list_fixedNegVSposRratio_keepCombination=list()
    for key in sorted(list(test_Phaly_posPPInum_fixedNegVSposRratio_keepCombination.keys())):
        #key=key/10
        print(key)
        test_Phaly_PPInum_list_fixedNegVSposRratio_keepCombination.append([key,test_Phaly_posPPInum_fixedNegVSposRratio_keepCombination[key],"P"])
        test_Phaly_PPInum_list_fixedNegVSposRratio_keepCombination.append([key,test_Phaly_negPPInum_fixedNegVSposRratio_keepCombination[key],"N"])
    test_Phaly_PPInum_frame_fixedNegVSposRratio_keepCombination=pd.DataFrame(test_Phaly_PPInum_list_fixedNegVSposRratio_keepCombination,columns=["-1Ratio","number","benchmark_status"])
    
    if ifplot:
        import matplotlib.pyplot as plt 
        import seaborn as sns
        plt.figure()
        ax = sns.barplot(x="-1Ratio", y="number",hue="benchmark_status", data=beforeFixRatio_test_Phaly_PPInum_frame)
        plt.show()
        plt.figure()
        ax = sns.barplot(x="-1Ratio", y="number",hue="benchmark_status", data=test_Phaly_PPInum_frame_fixedNegVSposRratio_keepCombination)
        plt.show()
    
    
    return(fixedNegVSposRratio_keepCombination_allIdx,
           beforeFixRatio_test_Phaly_PPInum_frame,
           test_Phaly_PPInum_frame_fixedNegVSposRratio_keepCombination)
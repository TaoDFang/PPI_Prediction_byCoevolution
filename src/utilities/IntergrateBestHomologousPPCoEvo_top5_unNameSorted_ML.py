import pandas as pd
import numpy as np
import copy
import os
import sys
import pickle
import multiprocessing as mp 


from collections import defaultdict
import matplotlib.pyplot as plt 


from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

#from xgboost import XGBClassifier

from sklearn.impute import SimpleImputer
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.impute import KNNImputer




# sys.path.append('/mnt/mnemo5/tao/code/MNF/src/tao_utilities/')

# from collect_topCoEvos import get_topRankingBetValue_dict
# from collect_topCoEvos import get_topRankingBetValue_dict_PosInSingleMSA
# from collect_topCoEvos import collect_topCoEvos_OnlyTopPosNeg
from collect_topCoEvos import get_topRanking_CoEvo_file


from EggNOGGroup_RelatedFuncs import return_cogPairs_fromproPair
from EggNOGGroup_RelatedFuncs import sepCogPairs_train_test_split

from biasCheck import phylaIntegration_replacingDCAScores
from biasCheck import phylaIntegration_replacingSubjectSpeDCAScores
from biasCheck import phylaIntegration_replacingOtherPhalaDCAScores

# from biasCheck import sort_arrayOtherPhyla
# from biasCheck import sort_arrayAllPhyla

from IntergrateBestHomologousPPCoEvo_unNameSorted import getMetaFrame_FullBestHomologousPP
from IntergrateBestHomologousPPCoEvo_unNameSorted import collect_BestHomologousDCAs_OneSpeOneScore_OnlyTopPosNeg

from read_benchmark_tableformatfile import getMetaFrame

from ML_parameters_setting import LR_withGridSearchCV_GroupKFold
from ML_parameters_setting import RF_withGridSearchCV_GroupKFold


def VariousReplacing_sepCogPairs_FullBestHomologousPP_BestHomologousDCAs_top5DCAs_uniquePhyla_ML_predictions(EggNOG_maxLevel,currentSpe_TaxID,STRING_Version,
                  Query_BestHomologousDCAs_dict,
                    EggNOG_group_level2,
                   ML_methods=["LR","RF"],
                  CoEvo_type="DCA",
                    topDCA_num=5,
                  deleting_column=None, # notice 0 colum are Query speceis 
                  DCA_thres=1,
                fillDCAValue=1,
                fillMissingValue=-1,
                replacing=True,
                replacing_strageties="",
                removeSamples_withManyNans=False,
                given_benchmark_folder=None,
                splitPosandNeg=True,
                sort_frame=True,
                CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING11.5/",
                prefix="",
                benchmark_suffix="STRINPhyPPI_Benchmark/",
                n_jobs=20,
                                                                                                    scoring_metrics=None,
                BestHomologousDCAs_dict_withStatus=True,):
    
    
    '''
    return prediciotn results of different machine leanning models;
    
    Here EggNOG_maxLevel,currentSpe_TaxID,STRING_Version acutally means Query id and Query max eggnog levels; 
    
    :param DCA_num :  number of top ranking DCA scores to be computated in data preparation step 
    :type DCA_num:  int 
    
    '''
    CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING"+STRING_Version+"/"
    input_root_folder=CoEvo_data_folder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_eggNOGfilteredData/"
    
    Subect_prefix="BestHomologousPPFor"+currentSpe_TaxID+"AtEggNOGmaxLevel"+EggNOG_maxLevel+"_"
    Benchmark_folder=input_root_folder+Subect_prefix+benchmark_suffix #

    
    ML_inputs=Benchmark_folder+"ML_inputs/"

    if not os.path.exists(ML_inputs):
        os.makedirs(ML_inputs)
        
        
    BestHomologousDCAs_predicted_results=dict()
        
    #here something is wrong ... whats problem. later remember to write down more details if we found problems 
    allPPI_allInfo_frame=getMetaFrame(EggNOG_maxLevel=EggNOG_maxLevel,currentSpe_TaxID=currentSpe_TaxID,
                                                      STRING_Version=STRING_Version,DCA_thres=DCA_thres,
                                                  given_benchmark_folder=given_benchmark_folder,
                                                  benchmark_suffix=benchmark_suffix,
                                                    splitPosandNeg=splitPosandNeg,
                                                  sort_frame=sort_frame,
                                                  CoEvo_data_folder=CoEvo_data_folder,
                                                 )
    
    allPPI_cogs=return_cogPairs_fromproPair(EggNOG_group_level2,currentSpe_TaxID,allPPI_allInfo_frame)
    
    
    OnlyTopPosNeg_NonPara_XBestHomologousDCAs,OnlyTopPosNeg_NonPara_YBestHomologousDCAs=collect_BestHomologousDCAs_OneSpeOneScore_OnlyTopPosNeg(allPPI_allInfo_frame,
                                                          Query_BestHomologousDCAs_dict,
                                                          ML_inputs,
                                                          CoEvo_type=CoEvo_type,
                                                        BestHomologousDCAs_dict_withStatus=BestHomologousDCAs_dict_withStatus)
    print("OnlyTopPosNeg_NonPara_XBestHomologousDCAs.shape:",OnlyTopPosNeg_NonPara_XBestHomologousDCAs.shape)
    #
    if removeSamples_withManyNans == True:
        NotManyNans_idx=np.sum(np.isnan(OnlyTopPosNeg_NonPara_XBestHomologousDCAs),axis=1)<=(len(list(Query_BestHomologousDCAs_dict.items())[0][1])-1-3)
        print(NotManyNans_idx)
        allPPI_cogs=[cogs for i,cogs  in enumerate(allPPI_cogs) if NotManyNans_idx[i]==True]
        allPPI_allInfo_frame=allPPI_allInfo_frame.loc[NotManyNans_idx,:]
        OnlyTopPosNeg_NonPara_XBestHomologousDCAs=OnlyTopPosNeg_NonPara_XBestHomologousDCAs[NotManyNans_idx,:]
        OnlyTopPosNeg_NonPara_YBestHomologousDCAs=OnlyTopPosNeg_NonPara_YBestHomologousDCAs[NotManyNans_idx]
        
        print("after removeSamples_withManyNans, OnlyTopPosNeg_NonPara_XBestHomologousDCAs.shape:",OnlyTopPosNeg_NonPara_XBestHomologousDCAs.shape)
        print("after removeSamples_withManyNans,len(allPPI_cogs):",len(allPPI_cogs))
    
    
    if deleting_column is not None:
        OnlyTopPosNeg_NonPara_XBestHomologousDCAs=np.delete(OnlyTopPosNeg_NonPara_XBestHomologousDCAs,deleting_column,axis=1)
    
    print("after deleting colum, OnlyTopPosNeg_NonPara_XBestHomologousDCAs.shape:",OnlyTopPosNeg_NonPara_XBestHomologousDCAs.shape)
    
  
    OnlyTopPosNeg_NonPara_XBestHomologousDCAs[np.isnan(OnlyTopPosNeg_NonPara_XBestHomologousDCAs)] = fillMissingValue
    if replacing == True:
        if replacing_strageties=="replacingDCAScores":
            print("The replacing_strageties is ", replacing_strageties)
            OnlyTopPosNeg_NonPara_XBestHomologousDCAs=phylaIntegration_replacingDCAScores(OnlyTopPosNeg_NonPara_XBestHomologousDCAs,
                                           fillDCAValue=fillDCAValue,
                                           fillMissingValue=fillMissingValue)

        elif replacing_strageties=="replacingSubjectSpeDCAScores":
            print("The replacing_strageties is ", replacing_strageties)
            OnlyTopPosNeg_NonPara_XBestHomologousDCAs=phylaIntegration_replacingSubjectSpeDCAScores(OnlyTopPosNeg_NonPara_XBestHomologousDCAs,
                                           fillDCAValue=fillDCAValue,
                                           fillMissingValue=fillMissingValue,
                                           SubjectFeaNum=topDCA_num,
                                            )

        elif replacing_strageties=="replacingOtherPhalaDCAScores":
            print("The replacing_strageties is ", replacing_strageties)
            OnlyTopPosNeg_NonPara_XBestHomologousDCAs=phylaIntegration_replacingOtherPhalaDCAScores(OnlyTopPosNeg_NonPara_XBestHomologousDCAs,
                                           fillDCAValue=fillDCAValue,
                                           fillMissingValue=fillMissingValue,
                                           SubjectFeaNum=topDCA_num,
                                            )
    
    
    #return(OnlyTopPosNeg_NonPara_XBestHomologousDCAs)

    BestHomologousDCAs_predicted_results["XBestHomologousDCAs"]=OnlyTopPosNeg_NonPara_XBestHomologousDCAs
    BestHomologousDCAs_predicted_results["YBestHomologousDCAs"]=OnlyTopPosNeg_NonPara_YBestHomologousDCAs
    
    cogs_train_array,XBestHomologousDCAs_train, XBestHomologousDCAs_test, yBestHomologousDCAs_train, yBestHomologousDCAs_test = sepCogPairs_train_test_split(allPPI_cogs,OnlyTopPosNeg_NonPara_XBestHomologousDCAs, OnlyTopPosNeg_NonPara_YBestHomologousDCAs, test_size=0.2, random_state=0)
    
    print("XBestHomologousDCAs_train.shape,yBestHomologousDCAs_train.shape,sum(yBestHomologousDCAs_train),yBestHomologousDCAs_test.shape,sum(yBestHomologousDCAs_test):",XBestHomologousDCAs_train.shape,yBestHomologousDCAs_train.shape,sum(yBestHomologousDCAs_train),yBestHomologousDCAs_test.shape,sum(yBestHomologousDCAs_test))
    
    
    
    BestHomologousDCAs_predicted_results["XBestHomologousDCAs_train"]=XBestHomologousDCAs_train
    BestHomologousDCAs_predicted_results["XBestHomologousDCAs_test"]=XBestHomologousDCAs_test
    BestHomologousDCAs_predicted_results["yBestHomologousDCAs_train"]=yBestHomologousDCAs_train
    BestHomologousDCAs_predicted_results["yBestHomologousDCAs_test"]=yBestHomologousDCAs_test
    
    
    if "LR" in ML_methods:
        print("train LR model now : ")

        clrBestHomologousDCAs=LR_withGridSearchCV_GroupKFold(XBestHomologousDCAs_train,yBestHomologousDCAs_train,cogs_train_array,n_jobs=n_jobs,scoring_metrics=scoring_metrics,)

        LR_yBestHomologousDCAs_predict=clrBestHomologousDCAs.predict(XBestHomologousDCAs_test)
        LR_yBestHomologousDCAs_predict_prob=clrBestHomologousDCAs.predict_proba(XBestHomologousDCAs_test)
        
        BestHomologousDCAs_predicted_results["LR"]=dict()
        BestHomologousDCAs_predicted_results["LR"]["Model"]=clrBestHomologousDCAs
        BestHomologousDCAs_predicted_results["LR"]["LR_yBestHomologousDCAs_predict"]=LR_yBestHomologousDCAs_predict
        BestHomologousDCAs_predicted_results["LR"]["LR_yBestHomologousDCAs_predict_prob"]=LR_yBestHomologousDCAs_predict_prob
        allPPI_allInfo_frame["LR_onesProb"]=clrBestHomologousDCAs.predict_proba(OnlyTopPosNeg_NonPara_XBestHomologousDCAs)[:,1]


    if "RF" in ML_methods:
        print("train RF model now : ")
        RFBestHomologousDCAs_model=RF_withGridSearchCV_GroupKFold(XBestHomologousDCAs_train,yBestHomologousDCAs_train,cogs_train_array,n_jobs=n_jobs,scoring_metrics=scoring_metrics,)
        RF_yBestHomologousDCAs_predict=RFBestHomologousDCAs_model.predict(XBestHomologousDCAs_test)
        RF_yBestHomologousDCAs_predict_prob=RFBestHomologousDCAs_model.predict_proba(XBestHomologousDCAs_test)
        
        BestHomologousDCAs_predicted_results["RF"]=dict()
        BestHomologousDCAs_predicted_results["RF"]["Model"]=RFBestHomologousDCAs_model
        BestHomologousDCAs_predicted_results["RF"]["RF_yBestHomologousDCAs_predict"]=RF_yBestHomologousDCAs_predict
        BestHomologousDCAs_predicted_results["RF"]["RF_yBestHomologousDCAs_predict_prob"]=RF_yBestHomologousDCAs_predict_prob
        allPPI_allInfo_frame["RF_onesProb"]=RFBestHomologousDCAs_model.predict_proba(OnlyTopPosNeg_NonPara_XBestHomologousDCAs)[:,1]
        
    BestHomologousDCAs_predicted_results["updated_allPPI_allInfo_frame"]=allPPI_allInfo_frame
        
        
    return BestHomologousDCAs_predicted_results
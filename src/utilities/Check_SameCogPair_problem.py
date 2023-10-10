import sys
import os 
from pathlib import Path
import pandas as pd
import multiprocessing as mp 
import numpy as np 

from collections import defaultdict
import copy

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier



sys.path.append('../src/utilities/')
from read_benchmark_tableformatfile import getMetaFrame


from EggNOGGroup_RelatedFuncs import return_cogPairs_fromproPair
from EggNOGGroup_RelatedFuncs import sepCogPairs_train_test_split


from collect_topCoEvos import get_topRanking_CoEvo_file
from collect_topCoEvos import collect_topCoEvos_OnlyTopPosNeg




from ML_parameters_setting import LR_withGridSearchCV
from ML_parameters_setting import RF_withGridSearchCV
from ML_parameters_setting import LR_withGridSearchCV_GroupKFold
from ML_parameters_setting import RF_withGridSearchCV_GroupKFold



def sepCogPairs_topDCAs_ML_predictions_getTestPPTuples(EggNOG_maxLevel,currentSpe_TaxID,STRING_Version,
                                       EggNOG_group_level2,
                  DCA_thres=1,filterIdx=None,
                    CoEvo_data_folder="CoEvo_data_STRING11.5/",
                    benchmark_suffix="STRINPhyPPI_Benchmark/"):
    '''

    
    '''
    input_root_folder=CoEvo_data_folder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_eggNOGfilteredData/"
    Benchmark_folder=input_root_folder+benchmark_suffix
    DCA_coevolutoin_path=input_root_folder+"coevolutoin_result_DCA/"
        
    topDCAs_predicted_results=dict()
        
        


    allPPI_allInfo_frame=getMetaFrame_withHighDCA(EggNOG_maxLevel=EggNOG_maxLevel,currentSpe_TaxID=currentSpe_TaxID,
                                                      STRING_Version="11.5",DCA_thres=DCA_thres,benchmark_suffix=benchmark_suffix)
    
    if filterIdx is not None:
        allPPI_allInfo_frame=allPPI_allInfo_frame.iloc[filterIdx,:]
    print("after filtering, ",allPPI_allInfo_frame.shape)
    

    allPPI_cogs=return_cogPairs_fromproPair(EggNOG_group_level2,currentSpe_TaxID,allPPI_allInfo_frame)

    
    allPPI_PPs=allPPI_allInfo_frame.loc[:,["STRING_ID1","STRING_ID2"]].values.tolist()
    allPPI_PPs=[tuple(pp) for pp in allPPI_PPs]
    print("len(allPPI_PPs):",len(allPPI_PPs))


  
    all_PPTuples_train,all_PPTuples_test = sepCogPairs_train_test_split_getTrainandTestPPTuples(allPPI_cogs,
                                                                                                                allPPI_PPs,
                                                                                                                test_size=0.20, 
                                                                                                                random_state=0)
    
    return(all_PPTuples_train,all_PPTuples_test)

def sepCogPairs_topDCAs_ML_predictions(EggNOG_maxLevel,currentSpe_TaxID,STRING_Version,
                                       EggNOG_group_level2,
                   ML_methods=["LR","RF"],
                    coevo_suffix="_pydcaFNAPC_array",
                  DCA_thres=1,DCA_number=50,selDca_number=20,
CoEvo_data_folder="CoEvo_data_STRING11.5/",
                                       prefix="",
                   benchmark_suffix="STRINPhyPPI_Benchmark/",
                                       downsample_prefix="",
                                       ifReCollect=False,
                                      saveFrame=False,
                                      overwrite=True,                                                                                      splitPosandNeg=True,
                                       sort_frame=True,
                                       use_multiprocessing=30,
                                        n_jobs=20,scoring_metrics=None,
                                      ):
    '''
    return prediciotn results of different machine leanning models;
    
    :param DCA_thres :  number of top ranking DCA scores to be computated in data preparation step 
    :type DCA_number:  int 
    
    :param DCA_number:  number of top ranking DCA scors used as input to Machine learning models
    :type selDca_number:  in t
    
    '''
    
    input_root_folder=CoEvo_data_folder+prefix+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_eggNOGfilteredData/"
    Benchmark_folder=input_root_folder+benchmark_suffix
        
    if coevo_suffix=="_pydcaFNAPC_array":
        DCA_coevolutoin_path=input_root_folder+downsample_prefix+"coevolutoin_result_DCA/"
        
        
    if coevo_suffix=="_alphafold_prob12":
        colab_output_path="/mnt/mnemo6/tao/colabfold_PPI_Coevolution/CoEvo_data_STRING"+STRING_Version+"/"
        DCA_coevolutoin_path=os.path.join(colab_output_path,f"EggNog_{currentSpe_TaxID}_EggNOGmaxLeve{EggNOG_maxLevel}/")




    ML_inputs=Benchmark_folder+"ML_inputs/"
    print("ML_inputs:",ML_inputs)

    if not os.path.exists(ML_inputs):
        os.makedirs(ML_inputs)
        
        
    topDCAs_predicted_results=dict()
        
        


    allPPI_allInfo_frame=getMetaFrame(EggNOG_maxLevel=EggNOG_maxLevel,currentSpe_TaxID=currentSpe_TaxID,
                                                      STRING_Version="11.5",DCA_thres=DCA_thres,
                                                 given_benchmark_folder=Benchmark_folder, benchmark_suffix=benchmark_suffix,splitPosandNeg=splitPosandNeg,sort_frame=sort_frame,
                                                  CoEvo_data_folder=CoEvo_data_folder,)
    

    allPPI_cogs=return_cogPairs_fromproPair(EggNOG_group_level2,currentSpe_TaxID,allPPI_allInfo_frame)
    
    
    allPPI_info=allPPI_allInfo_frame.loc[:,["STRING_ID1","STRING_ID2","len1","len2"]].values.tolist()
    print("len(allPPI_info):",len(allPPI_info))

    
    
    #here use suffix to decide if collect dca or alphafold 
    
    if coevo_suffix=="_pydcaFNAPC_array":
        topRanking_pydcaFNAPC_file=Benchmark_folder+"DCA_thres_"+str(DCA_thres)+"_topRanking_pydcaFNAPC_frame.csv" 
        
    if coevo_suffix=="_alphafold_prob12":
        topRanking_pydcaFNAPC_file=Benchmark_folder+"DCA_thres_"+str(DCA_thres)+"_topRanking_alphafoldProb12_frame.csv" 
        
        
    if coevo_suffix=="_alphafold_inversecontactmap":
        topRanking_pydcaFNAPC_file=Benchmark_folder+"DCA_thres_"+str(DCA_thres)+"_topRanking_alphafoldInverseCm_frame.csv" 
        
    top_pydcaFNAPC_frame=get_topRanking_CoEvo_file(topRanking_CoEvo_file=topRanking_pydcaFNAPC_file,
                                                       coevolutoin_path=DCA_coevolutoin_path,
                                                       coevo_suffix=coevo_suffix,
                                                       allPPI_info=allPPI_info, 
                                                       returnDic=False,
                                                       overwrite=overwrite,
                                                  use_multiprocessing=use_multiprocessing)

    if coevo_suffix=="_pydcaFNAPC_array":
        OnlyTopPosNeg_NonPara_XtopDCAs,OnlyTopPosNeg_NonPara_YtopDCAs=collect_topCoEvos_OnlyTopPosNeg(top_pydcaFNAPC_frame,
         allPPI_allInfo_frame,"DCA"+"DCA_thres_"+str(DCA_thres),DCA_thres,DCA_number,selDca_number,ML_inputs,ifReCollect)
    else:
        OnlyTopPosNeg_NonPara_XtopDCAs,OnlyTopPosNeg_NonPara_YtopDCAs=collect_topCoEvos_OnlyTopPosNeg(top_pydcaFNAPC_frame,
            allPPI_allInfo_frame,coevo_suffix[1:]+"_",DCA_thres,
            DCA_number,selDca_number,ML_inputs,ifReCollect)
    
    topDCAs_predicted_results["XtopDCAs"]=OnlyTopPosNeg_NonPara_XtopDCAs
    topDCAs_predicted_results["YtopDCAs"]=OnlyTopPosNeg_NonPara_YtopDCAs
    
  
    cogs_train_array,XtopDCAs_train, XtopDCAs_test, ytopDCAs_train, ytopDCAs_test = sepCogPairs_train_test_split(allPPI_cogs,OnlyTopPosNeg_NonPara_XtopDCAs, OnlyTopPosNeg_NonPara_YtopDCAs, test_size=0.20, random_state=0)
 
    print("XtopDCAs_train.shape,ytopDCAs_train.shape,sum(ytopDCAs_train),ytopDCAs_test.shape,sum(ytopDCAs_test):",XtopDCAs_train.shape,ytopDCAs_train.shape,sum(ytopDCAs_train),ytopDCAs_test.shape,sum(ytopDCAs_test))
    
    
    topDCAs_predicted_results["XtopDCAs_train"]=XtopDCAs_train
    topDCAs_predicted_results["XtopDCAs_test"]=XtopDCAs_test
    topDCAs_predicted_results["ytopDCAs_train"]=ytopDCAs_train
    topDCAs_predicted_results["ytopDCAs_test"]=ytopDCAs_test
    
    
    
    
    if "LR" in ML_methods:
        print("train LR now _:")
        clrtopDCAs = LR_withGridSearchCV_GroupKFold(XtopDCAs_train,ytopDCAs_train,cogs_train_array,n_jobs=n_jobs,scoring_metrics=scoring_metrics)

        LR_ytopDCAs_predict=clrtopDCAs.predict(XtopDCAs_test)
        LR_ytopDCAs_predict_prob=clrtopDCAs.predict_proba(XtopDCAs_test)
        
        topDCAs_predicted_results["LR"]=dict()
        topDCAs_predicted_results["LR"]["Model"]=clrtopDCAs
        topDCAs_predicted_results["LR"]["LR_ytopDCAs_predict"]=LR_ytopDCAs_predict
        topDCAs_predicted_results["LR"]["LR_ytopDCAs_predict_prob"]=LR_ytopDCAs_predict_prob
        allPPI_allInfo_frame["LR_onesProb"]=clrtopDCAs.predict_proba(OnlyTopPosNeg_NonPara_XtopDCAs)[:,1]
        
        #topDCAs_predicted_results["FullDatasetPredictedProb"]["LR"]=clrtopDCAs.predict_proba(OnlyTopPosNeg_NonPara_XtopDCAs)

    if "RF" in ML_methods:
        print("train RF now _:")
        RFtopDCAs_model = RF_withGridSearchCV_GroupKFold(XtopDCAs_train,ytopDCAs_train,cogs_train_array,n_jobs=n_jobs,scoring_metrics=scoring_metrics)
        
        RF_ytopDCAs_predict=RFtopDCAs_model.predict(XtopDCAs_test)
        RF_ytopDCAs_predict_prob=RFtopDCAs_model.predict_proba(XtopDCAs_test)
        
        topDCAs_predicted_results["RF"]=dict()
        topDCAs_predicted_results["RF"]["Model"]=RFtopDCAs_model
        topDCAs_predicted_results["RF"]["RF_ytopDCAs_predict"]=RF_ytopDCAs_predict
        topDCAs_predicted_results["RF"]["RF_ytopDCAs_predict_prob"]=RF_ytopDCAs_predict_prob
        allPPI_allInfo_frame["RF_onesProb"]=RFtopDCAs_model.predict_proba(OnlyTopPosNeg_NonPara_XtopDCAs)[:,1]
        #topDCAs_predicted_results["FullDatasetPredictedProb"]["RF"]=RFtopDCAs_model.predict_proba(OnlyTopPosNeg_NonPara_XtopDCAs)
    topDCAs_predicted_results["updated_allPPI_allInfo_frame"]=allPPI_allInfo_frame
    
    save_allPPI_allInfo_file=ML_inputs+"DCA_thres_"+str(DCA_thres)+"sepCogPairs_FullDataset_ML_onesProb.csv"
    print("save_allPPI_allInfo_file:",save_allPPI_allInfo_file)
    if saveFrame:
        allPPI_allInfo_frame.to_csv(save_allPPI_allInfo_file,
                                    header=True,index=None,sep="\t")
    return topDCAs_predicted_results 
    
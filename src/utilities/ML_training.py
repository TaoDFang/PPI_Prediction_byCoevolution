import os 
import numpy as np
import pandas as pd
import sys



from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

from  EggNOGGroup_RelatedFuncs import return_cogPairs_fromproPair
from  EggNOGGroup_RelatedFuncs import sepCogPairs_train_test_split
from EggNOGGroup_RelatedFuncs import sepCogPairs_train_test_split_getTrainandTestPPTuples

from ML_parameters_setting import LR_withGridSearchCV_GroupKFold
from ML_parameters_setting import RF_withGridSearchCV_GroupKFold



def sepCogPairs_ML_predictions_allTypeFeas(currentSpe_TaxID,
                                       EggNOG_group_level2,
                allPPI_allInfo_frame,
                inputFea_matrixList,
                Y_matrix, 
                ML_methods=["LR","RF"],
                filterIdx=None,
                n_jobs=20,
                scoring_metrics=None,
                getTrainandTestPPTuples=False,):
    '''
    
    '''
    allPPI_allInfo_frame_deepcp=allPPI_allInfo_frame.copy(deep=True)
        
    topFeas_predicted_results=dict()
        

    if len(inputFea_matrixList)==1:
        X_matrix=inputFea_matrixList[0]
        if X_matrix.ndim==1:
            X_matrix=np.reshape(X_matrix,(-1,1))
    else:
    
        X_matrix=np.concatenate((inputFea_matrixList), axis=1)
        
        
    if filterIdx is not None:
        X_matrix=X_matrix[filterIdx,:]
        Y_matrix=Y_matrix[filterIdx]
        allPPI_allInfo_frame_deepcp=allPPI_allInfo_frame_deepcp.iloc[filterIdx,:]
    print("after filtering, ",X_matrix.shape,Y_matrix.shape)
    
    topFeas_predicted_results["XtopFeas"]=X_matrix
    topFeas_predicted_results["YtopFeas"]=Y_matrix
    
    allPPI_cogs=return_cogPairs_fromproPair(EggNOG_group_level2,currentSpe_TaxID,allPPI_allInfo_frame_deepcp)
    if getTrainandTestPPTuples==True:
        all_PPTuples=allPPI_allInfo_frame_deepcp.loc[:,["STRING_ID1","STRING_ID2"]].values.tolist()
        all_PPTuples=[tuple(pp) for pp in all_PPTuples]
        all_PPTuples_train,all_PPTuples_test=sepCogPairs_train_test_split_getTrainandTestPPTuples(allPPI_cogs,all_PPTuples,test_size=0.20, random_state=0)
        return(all_PPTuples_train,all_PPTuples_test)
    else:
    
        cogs_train_array,XtopFeas_train, XtopFeas_test, ytopFeas_train, ytopFeas_test = sepCogPairs_train_test_split(allPPI_cogs,
                                                                                                    X_matrix, 
                                                                                                    Y_matrix, 
                                                                                                    test_size=0.20, random_state=0)

        print("XtopFeas_train.shape,ytopFeas_train.shape,sum(ytopFeas_train),ytopFeas_test.shape,sum(ytopFeas_test):",XtopFeas_train.shape,ytopFeas_train.shape,sum(ytopFeas_train),ytopFeas_test.shape,sum(ytopFeas_test))


        topFeas_predicted_results["XtopFeas_train"]=XtopFeas_train
        topFeas_predicted_results["XtopFeas_test"]=XtopFeas_test
        topFeas_predicted_results["ytopFeas_train"]=ytopFeas_train
        topFeas_predicted_results["ytopFeas_test"]=ytopFeas_test
    
    
    
    
    if "LR" in ML_methods:
        print("train LR model now : ")
        clrtopFeas=LR_withGridSearchCV_GroupKFold(XtopFeas_train,ytopFeas_train,cogs_train_array,
                                                 n_jobs=n_jobs,scoring_metrics=scoring_metrics)

        LR_ytopFeas_predict=clrtopFeas.predict(XtopFeas_test)
        LR_ytopFeas_predict_prob=clrtopFeas.predict_proba(XtopFeas_test)
        
        topFeas_predicted_results["LR"]=dict()
        topFeas_predicted_results["LR"]["Model"]=clrtopFeas
        topFeas_predicted_results["LR"]["LR_ytopFeas_predict"]=LR_ytopFeas_predict
        topFeas_predicted_results["LR"]["LR_ytopFeas_predict_prob"]=LR_ytopFeas_predict_prob
        
        LR_ytopFeas_predict_train=clrtopFeas.predict(XtopFeas_train)
        LR_ytopFeas_predict_prob_train=clrtopFeas.predict_proba(XtopFeas_train)
        topFeas_predicted_results["LR"]["LR_ytopFeas_predict_train"]=LR_ytopFeas_predict_train
        topFeas_predicted_results["LR"]["LR_ytopFeas_predict_prob_train"]=LR_ytopFeas_predict_prob_train
        
        
        allPPI_allInfo_frame_deepcp["LR_onesProb"]=clrtopFeas.predict_proba(X_matrix)[:,1]

    if "RF" in ML_methods:
        print("train RF model now : ")
        RFtopFeas_model=RF_withGridSearchCV_GroupKFold(XtopFeas_train,ytopFeas_train,cogs_train_array,
                                                      n_jobs=n_jobs,scoring_metrics=scoring_metrics)
        RF_ytopFeas_predict=RFtopFeas_model.predict(XtopFeas_test)
        RF_ytopFeas_predict_prob=RFtopFeas_model.predict_proba(XtopFeas_test)
        
        topFeas_predicted_results["RF"]=dict()
        topFeas_predicted_results["RF"]["Model"]=RFtopFeas_model
        topFeas_predicted_results["RF"]["RF_ytopFeas_predict"]=RF_ytopFeas_predict
        topFeas_predicted_results["RF"]["RF_ytopFeas_predict_prob"]=RF_ytopFeas_predict_prob
        
        RF_ytopFeas_predict_train=RFtopFeas_model.predict_proba(XtopFeas_train)
        topFeas_predicted_results["RF"]["RF_ytopFeas_predict_train"]=RF_ytopFeas_predict_train
        RF_ytopFeas_predict_prob_train=RFtopFeas_model.predict_proba(XtopFeas_train)
        topFeas_predicted_results["RF"]["RF_ytopFeas_predict_prob_train"]=RF_ytopFeas_predict_prob_train
        
        allPPI_allInfo_frame_deepcp["RF_onesProb"]=RFtopFeas_model.predict_proba(X_matrix)[:,1]
    topFeas_predicted_results["updated_allPPI_allInfo_frame"]=allPPI_allInfo_frame_deepcp

    return topFeas_predicted_results
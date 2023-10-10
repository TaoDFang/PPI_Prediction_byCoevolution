import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import random 
import copy 

from typing import Optional

from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve






def F1_score(PPV,TPR):
    return(2*(PPV*TPR)/(PPV+TPR))

def accuracy(TP,TN,P,N):
    return((TP+TN)/(P+N))

def DCA_RocCurve(X_test,Y_test, count_label="rate",legend=None,show_legend=True,step=1000,zoom_thres=None,tiny_randomnoise_range=0,color=None,return_metrics=False,
    figure_handle=plt,ifplot=True,plotType="ROC",
                randomplot=False,
                fontdict=None):
    """
    Here draw ROC prediction performce plot ;

    :param X_test : predicted probabilities
    :type X_test: np.array either with one column or multiple columns but prediction probability have to be in last column 

    :param Y_test : binary true labels
    :type Y_test: np.array 

    :param count_label : "rate" or "count " ; either draw roc curve in rate or real count
    :type count_label: string  

    :param legend : text used to label prediction curve
    :type legend: string

    :param step : the nummber of  seletcted postive predictions in each step   from ranked predicted probabities; small value take more time to plot
    :type step: int 

    :param zoom_thres : use to zoom into begining part of ROC curver; from 0 to number of sampels 
    :type zoom_thres: int
    
    :param tiny_randomnoise_range : use to random sort same values 
    :type tiny_randomnoise_range: flot

    :return: plot a ROC curve
    :rtype: None

    """
    
    TPR_list=list()
    FPR_list=list()
    
    
    X_test=np.array(X_test)
    Y_test=np.array(Y_test)
    
    if X_test.ndim>1:
        X_DCA=X_test[:,-1]
    else:
        X_DCA=X_test
        
    P=sum(Y_test==1)
    N=sum(Y_test==0)
    
    if not randomplot:  # only draw random plot withouth zoom in         
        np.random.seed(0)
        noise=np.random.normal(0,tiny_randomnoise_range,X_DCA.shape[0])
        X_DCA=noise+X_DCA

        Ascending_orderIdx=np.argsort(X_DCA)
        Descending_orderIdx=Ascending_orderIdx[::-1]
        X_DCA_ori=X_DCA[Descending_orderIdx]
        Y_test_ori=Y_test[Descending_orderIdx]

        if zoom_thres is not None:
            #X_DCA=X_DCA_ori[0:zoom_thres]
            Y_test=Y_test_ori[0:zoom_thres]
        else:
            #X_DCA=X_DCA_ori
            Y_test=Y_test_ori
    else: 
        Y_test=copy.deepcopy(Y_test)
        random.Random(10).shuffle(Y_test)


    
    TPR_step_list=[0]
    FPR_step_list=[0]
    for i in range(step,len(X_DCA),step):
        y_top_list_rank=Y_test[(i-step):(i)]

        TPR_step=sum([1  for  r in y_top_list_rank if  r==1])
        FPR_step=sum([1  for  r in y_top_list_rank if  r==0])
        TPR=sum(TPR_step_list)+TPR_step
        FPR=sum(FPR_step_list)+FPR_step
        
        TPR_step_list.append(TPR_step)
        FPR_step_list.append(FPR_step)
  
        TPR_list.append(TPR)
        FPR_list.append(FPR) 
    
    TPR_rate_list =[TP/P for TP in TPR_list]
    FPR_rate_list =[FP/N for FP in FPR_list]
        
    PPV_rate_list= [TPR_list[i]/(TPR_list[i]+FPR_list[i]) for i in range(len(TPR_list))] 
        
    if ifplot:
        if plotType=="ROC":
            if count_label=="count":
                figure_handle.plot(FPR_list,TPR_list,label=legend,color=color)
                try: # for plt
                    figure_handle.xlabel("False positive count",fontdict=fontdict)
                    figure_handle.ylabel("True positive count",fontdict=fontdict)
                except: # for plt.subplot
                    figure_handle.set_xlabel("False positive count",fontdict=fontdict)
                    figure_handle.set_ylabel("True positive count",fontdict=fontdict)
                    
            elif count_label=="rate":
                figure_handle.plot(FPR_rate_list,TPR_rate_list,label=legend,color=color)
                try:
                    figure_handle.xlabel("False positive rate",fontdict=fontdict)
                    figure_handle.ylabel("True positive rate",fontdict=fontdict)
                except:
                    figure_handle.set_xlabel("False positive rate",fontdict=fontdict)
                    figure_handle.set_ylabel("True positive rate",fontdict=fontdict)
                
        elif plotType=="PR":
            if count_label=="count":
                pass
            elif count_label=="rate":
                figure_handle.plot(TPR_rate_list,PPV_rate_list,label=legend,color=color)
                try:
                    figure_handle.xlabel("Recall",fontdict=fontdict)
                    figure_handle.ylabel("Precision",fontdict=fontdict)
                except:
                    figure_handle.set_xlabel("Recall",fontdict=fontdict)
                    figure_handle.set_ylabel("Precision",fontdict=fontdict)

                

        if show_legend:
            figure_handle.legend(loc="lower right")
    
        
    if return_metrics:
        
        ROC_AUC=np.round(np.trapz(TPR_rate_list, x=FPR_rate_list),4)
        
        PR_AUC=np.round(np.trapz(PPV_rate_list, x=TPR_rate_list),4)
        

        F1_score_lastStep=F1_score(PPV_rate_list[-1],TPR_rate_list[-1])


        

        TNR_list=[N-FP for FP in FPR_list]
        ACC_lastStep=accuracy(TPR_list[-1],TNR_list[-1],P,N)
        return([ROC_AUC,PR_AUC,F1_score_lastStep,ACC_lastStep])


def Random_RocCurve(ori_Y_test, count_label="rate",legend=None,show_legend=True,step=1000,zoom_thres=None,
                   color=None):
    """
    Here draw ROC prediction performce plot ;

    :param X_test : predicted probabilities
    :type X_test: np.array either with one column or multiple columns but prediction probability have to be in last column 

    :param Y_test : binary true labels
    :type Y_test: np.array 

    :param count_label : "rate" or "count " ; either draw roc curve in rate or real count
    :type count_label: string  

    :param legend : text used to label prediction curve
    :type legend: string

    :param step : the nummber of  seletcted postive predictions in each step   from ranked predicted probabities; small value take more time to plot
    :type step: int 

    :param zoom_thres : use to zoom into begining part of ROC curver; from 0 to number of sampels 
    :type zoom_thres: int

    :return: plot a ROC curve
    :rtype: None

    """
    
    TPR_list=list()
    FPR_list=list()
    
    
    Y_test=copy.deepcopy(ori_Y_test)
    random.Random(10).shuffle(Y_test)

    
    if zoom_thres is not None:
        Y_test=Y_test[0:zoom_thres]
    

    TPR_step_list=[0]
    FPR_step_list=[0]
    for i in range(step,len(Y_test),step):
        y_top_list_rank=Y_test[(i-step):(i)]

        TPR_step=sum([1  for  r in y_top_list_rank if  r==1])
        FPR_step=sum([1  for  r in y_top_list_rank if  r==0])
        TPR=sum(TPR_step_list)+TPR_step
        FPR=sum(FPR_step_list)+FPR_step
        
        TPR_step_list.append(TPR_step)
        FPR_step_list.append(FPR_step)
  
        TPR_list.append(TPR)
        FPR_list.append(FPR) 
    
    P=sum(Y_test==1)
    N=sum(Y_test==0)
    TPR_rate_list =[TP/P for TP in TPR_list]
    FPR_rate_list =[FP/N for FP in FPR_list]
        


    if count_label=="count":
        plt.plot(FPR_list,TPR_list,label=legend,color=color)
    elif count_label=="rate":
        plt.plot(FPR_rate_list,TPR_rate_list,label=legend,color=color)
        
    if show_legend:
        plt.legend(loc="lower right")
    plt.xlabel("FPR")
    plt.ylabel("TRP")

    
    
    
def get_STRINGScoreFromDict(inputPPList,allSTRING_dict):
    
    returned_score=[allSTRING_dict[(spp0,spp1)] if \
                                (spp0,spp1) in allSTRING_dict\
                                   else 0
                                for spp0,spp1 in inputPPList ]
    returned_score=np.array(returned_score)
    returned_score.reshape(-1, 1) # 1 column 
    return(returned_score)
    
    
    

    
def get_F1_list(PPV_list,TPR_list):
    F1_list=list()
    for i in range(len(PPV_list)):
        PPV,TPR=PPV_list[i],TPR_list[i]
        if (PPV+TPR)!=0:
            F1=2*(PPV*TPR)/(PPV+TPR)
        else:
            F1=0
        F1_list.append(F1)
    return(F1_list)

def show_bestF1_fromPrecisionAndRecall(PPV_list,TPR_list,P,return_bestF1=False):
    F1_list=get_F1_list(PPV_list,TPR_list)
        
    F1_max=max(F1_list)
    F1_max_idx=F1_list.index(F1_max)
    BestF1,BestPrecision,BestRecall,BestTP=F1_list[F1_max_idx],PPV_list[F1_max_idx],TPR_list[F1_max_idx],TPR_list[F1_max_idx]*P
    print(f"When achieve max F1, precision, recall,TP: ",
          np.round(BestF1,decimals=4),np.round(BestPrecision,decimals=4),np.round(BestRecall,decimals=4),np.round(BestTP,decimals=4),)
    if return_bestF1:
        return(BestF1,BestPrecision,BestRecall,BestTP)
        
    

    
def Precison_Recall_RocCurve_OnIndepedentBenchmark(Y_Prediction,Y_Prediction_PPs,IndepedentBenchmark_pps_dict, count_label="rate",legend=None,step=1000,zoom_thres=None,tiny_randomnoise_range=0,print_results=False):
    # here pps in IndepedentBenchmark_pps_dict are consider actuall positive 
    
    TPR_list=list() # recall 
    PPV_list=list() # precision 
    #ident = [0.0, 1.0]
    Y_Prediction=np.array(Y_Prediction)
    Y_Prediction_PPs=np.array(Y_Prediction_PPs)

    np.random.seed(0)
    noise=np.random.normal(0,tiny_randomnoise_range,Y_Prediction.shape[0])
    Y_Prediction=noise+Y_Prediction
    
    Ascending_orderIdx=np.argsort(Y_Prediction)
    Descending_orderIdx=Ascending_orderIdx[::-1]
    
    Y_Prediction=Y_Prediction[Descending_orderIdx]
    Y_Prediction_PPs=Y_Prediction_PPs[Descending_orderIdx]
    
    if zoom_thres is not None:
        Y_Prediction=Y_Prediction[0:zoom_thres]
        Y_Prediction_PPs=Y_Prediction_PPs[0:zoom_thres]
    
    P=len(IndepedentBenchmark_pps_dict)
    

    TPR_step_list=[0]
    for i in range(step,len(Y_Prediction_PPs),step):
        Y_Prediction_PPs_rank=Y_Prediction_PPs[(i-step):(i)]

        #TPR=0
        TPR_step=sum([1  for  p1,p2 in Y_Prediction_PPs_rank if ((p1, p2) in IndepedentBenchmark_pps_dict) or ((p2, p1) in IndepedentBenchmark_pps_dict) ])
        TPR=sum(TPR_step_list)+TPR_step
        
        TPR_step_list.append(TPR_step)
        
        if count_label=="rate":
            TPR,PPV=TPR/P,TPR/i
        elif count_label=="count":
            TPR,PPV=TPR,TPR/i
        #print(TPR)
        TPR_list.append(TPR)
        PPV_list.append(PPV)  
        


    plt.plot(TPR_list,PPV_list,label=legend)
    plt.legend(loc="lower right")
    plt.xlabel("Recall/TPR")
    plt.ylabel("Precision/PPV")
    
    if print_results:
        return(PPV_list,TPR_list)
        
def Precison_Recall_RocCurve_OnlyGetRecallAndPrecision_OnIndepedentBenchmark(Pos_Prediction,Pos_Prediction_PPs,IndepedentBenchmark_pps_dict, count_label="rate",legend=None,step=1000,zoom_thres=None, filter_thres=None,
                                                                             tiny_randomnoise_range=0,print_results=False,return_bestF1=False):
    # Here Pos_Prediction only contain positive predictions , thats top predictions 
    TPR_list=list() # recall 
    PPV_list=list() # precision 

    Pos_Prediction=np.array(Pos_Prediction)
    Pos_Prediction_PPs=np.array(Pos_Prediction_PPs)

    np.random.seed(0)
    noise=np.random.normal(0,tiny_randomnoise_range,Pos_Prediction.shape[0])
    Pos_Prediction=noise+Pos_Prediction
    
    Ascending_orderIdx=np.argsort(Pos_Prediction)
    Descending_orderIdx=Ascending_orderIdx[::-1]
    
    Pos_Prediction=Pos_Prediction[Descending_orderIdx]
    Pos_Prediction_PPs=Pos_Prediction_PPs[Descending_orderIdx]
    
    if zoom_thres is not None:
        Pos_Prediction=Pos_Prediction[0:zoom_thres]
        Pos_Prediction_PPs=Pos_Prediction_PPs[0:zoom_thres]
        
    if filter_thres is not None:
        Pos_Prediction_old=copy.deepcopy(Pos_Prediction)
        Pos_Prediction=Pos_Prediction[Pos_Prediction_old>filter_thres]
        Pos_Prediction_PPs=Pos_Prediction_PPs[Pos_Prediction_old>filter_thres]
    
    P=len(IndepedentBenchmark_pps_dict)
     
    TPR_step_list=[0]
    for i in range(step,len(Pos_Prediction_PPs),step):
        Pos_Prediction_PPs_rank=Pos_Prediction_PPs[(i-step):(i)]

        #TPR=0
        TPR_step=sum([1  for  p1,p2 in Pos_Prediction_PPs_rank if ((p1, p2) in IndepedentBenchmark_pps_dict) or ((p2, p1) in IndepedentBenchmark_pps_dict) ])
        TPR=sum(TPR_step_list)+TPR_step
        
        TPR_step_list.append(TPR_step)
        
        if count_label=="rate":
            TPR,PPV=TPR/P,TPR/i
        elif count_label=="count":
            pass
        TPR_list.append(TPR)
        PPV_list.append(PPV)  

    plt.plot(TPR_list,PPV_list,label=legend)
    plt.legend(loc="lower right")
    plt.xlabel("Recall/TPR")
    plt.ylabel("Precision/PPV")
    
    if print_results:
        print(f"{legend}:postiveControlCount:{P};allPredCount:{len(Pos_Prediction_PPs)};")
        
        F1_list=get_F1_list(PPV_list,TPR_list)
        FinalF1,FinalPrecision,FinalRecall,FinalTP=F1_list[-1],PPV_list[-1],TPR_list[-1],TPR_list[-1]*P
        print(f"Final F1, precision, recall,TP: ",
        np.round(FinalF1,decimals=4),np.round(FinalPrecision,decimals=4),np.round(FinalRecall,decimals=4),np.round(FinalTP,decimals=4),)
        if return_bestF1:
            BestF1,BestPrecision,BestRecall,BestTP=show_bestF1_fromPrecisionAndRecall(PPV_list,TPR_list,P,return_bestF1)
            return(FinalF1,FinalPrecision,FinalRecall,FinalTP,BestF1,BestPrecision,BestRecall,BestTP)
        
    
def countPP_overlap(input_PPs,benchmarkList):
    count=0
    for p1,p2 in input_PPs:
        if ((p1,p2) in benchmarkList) or ((p2,p1) in benchmarkList):
            count+=1
    print(count," / ", len(benchmarkList),":",count/len(benchmarkList))
    
    
def DCA_RocCurve_fromDataframe(DataFrame:pd.DataFrame,
                               y_label:str,
                               x_label_list:list,
                               count_label: str ="rate",
                                step:int=10,
                                plotType:str="ROC",
                               legend_suffix:Optional[str] = None,
                              ) -> None:
    # https://docs.python.org/3/library/typing.html
    
    for x_label in x_label_list:
        DCA_RocCurve(np.array(DataFrame[x_label]), 
                     np.array(DataFrame[y_label]),
                     count_label=count_label,legend=f"{x_label};{legend_suffix}",step=step,
                    plotType=plotType)
    plt.show()
    
    
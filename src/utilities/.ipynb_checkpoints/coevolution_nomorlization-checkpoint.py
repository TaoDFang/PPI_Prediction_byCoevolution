import numpy as np


def get_idx2protein_dict(pp_dict):
    Unique_Pro1=set([p1 for p1 , p2 in pp_dict.keys()])
    print(len(Unique_Pro1))
    Unique_Pro2=set([p2 for p1 , p2 in pp_dict.keys()])
    print(len(Unique_Pro2))
    Unique_Pro=sorted(list(Unique_Pro1.union(Unique_Pro2)))
    print(len(Unique_Pro))
    print(Unique_Pro[0:3])
    Unique_idx2pro_dict={i:p for i,p in enumerate(Unique_Pro)}
    return(Unique_Pro,Unique_idx2pro_dict)
    
    
def get_proteinLevel_pairwise_matrix_dig(Unique_Pro,Unique_idx2pro_dict,top_coevo_dict,topDCA_num,):
    topDCA_idx=[3*i for i in range(topDCA_num)]
    pairwise_matrix_Dig=np.empty((len(Unique_idx2pro_dict)*topDCA_num,len(Unique_idx2pro_dict)*topDCA_num))
    pairwise_matrix_Dig[:] = np.nan
    print(pairwise_matrix_Dig.shape)
    
    count=0
    for i in range(len(Unique_Pro)):
        for j in range(i+1,len(Unique_Pro)):
    #     for j in range(len(Unique_Ecoli_Pro)): # this one does not work because our pp record is sorted 
            pro_i,pro_j=Unique_idx2pro_dict[i],Unique_idx2pro_dict[j]
            if (pro_i,pro_j) in top_coevo_dict:
                count +=1
                pairwise_matrix_Dig[i*topDCA_num,(topDCA_num*j):(topDCA_num*j+topDCA_num)]=[top_coevo_dict[(pro_i,pro_j)][ii] for ii in topDCA_idx]
                #pairwise_matrix_Dig[i,j]=top_coevo_dict[(pro_i,pro_j)][0]
    print(count)
    return(pairwise_matrix_Dig)

def get_proteinLevel_pairwise_matrix_full(Unique_Pro,Unique_idx2pro_dict,top_coevo_dict,topDCA_num,):
    pairwise_matrix_Dig=get_proteinLevel_pairwise_matrix_dig(Unique_Pro,Unique_idx2pro_dict,top_coevo_dict,topDCA_num,)

    pairwise_matrix= np.nan_to_num(pairwise_matrix_Dig.T) + np.nan_to_num(pairwise_matrix_Dig)
    pairwise_matri_mask=np.ma.masked_array(np.nan_to_num(pairwise_matrix), mask=(pairwise_matrix==0) )
    pairwise_matrix=pairwise_matri_mask.filled(np.nan)
    
    return pairwise_matrix



def proteinLevel_APC(pairwise_matrix):
    pairwise_matrix_sum=np.nansum(pairwise_matrix)
    pairwise_matrix_rowSum=np.nansum(pairwise_matrix,axis=1)
    pairwise_matrix_colSum=np.nansum(pairwise_matrix,axis=0)
    pairwise_matrix_rowCount=np.sum(~np.isnan(pairwise_matrix),axis=1)
    pairwise_matrix_colCount=np.sum(~np.isnan(pairwise_matrix),axis=0)
    pairwise_matrix_ave=pairwise_matrix_sum/(np.sum(pairwise_matrix_rowCount))

    print("pairwise_matrix_sum,pairwise_matrix_ave",pairwise_matrix_sum,pairwise_matrix_ave)
    print("sum(pairwise_matrix_rowCount)",sum(pairwise_matrix_rowCount))
    print("sum(pairwise_matrix_colCount)",sum(pairwise_matrix_colCount))
    print("pairwise_matrix.shape:",pairwise_matrix.shape)
    print("pairwise_matrix.shape[0]*pairwise_matrix.shape[1]:",pairwise_matrix.shape[0]*pairwise_matrix.shape[1])

    # add assert here ? it should never has 0s row/colum count
    print("np.where(pairwise_matrix_rowCount==0):",np.where(pairwise_matrix_rowCount==0))


    APCed_pairwise_matrix=np.empty((pairwise_matrix.shape[0],pairwise_matrix.shape[1]))
    APCed_pairwise_matrix[:] = np.nan

    for index in range(pairwise_matrix.shape[0]):
        for column in range(pairwise_matrix.shape[1]):
            APCed_pairwise_matrix[index,column]=pairwise_matrix[index,column]-\
            (pairwise_matrix_rowSum[index]/pairwise_matrix_rowCount[index])*(pairwise_matrix_colSum[column]/pairwise_matrix_colCount[column])/(pairwise_matrix_ave)
    
    return(APCed_pairwise_matrix)
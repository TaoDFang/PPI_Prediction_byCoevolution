from pydca.meanfield_dca import meanfield_dca # this need to load in py37_pydca enviroment
import numpy as np 



def post_pydcaFN2array_APC(QueryPro_pro1,QueryPro_pro2,L1,L2,DCA_ori,DCA_folder):
    DCA_array=np.zeros((L1+L2,L1+L2))

    array_rowIdx=np.array([t[0] for t in DCA_ori[:,0]])
    array_colIdx=np.array([t[1] for t in DCA_ori[:,0]])
    array_values=DCA_ori[:,1]
    DCA_array[array_rowIdx,array_colIdx]=array_values

    msa_len=L1+L2
    DCA_array_full = DCA_array.T + DCA_array
    np.fill_diagonal(DCA_array_full,DCA_array.diagonal())# as MI of column itself is not  0!!!!


    DCA_array_full_sum=DCA_array_full.sum()/ (msa_len*msa_len-msa_len)   #DCA_array_full.sum()
    DCA_array_full_rowSum=DCA_array_full.sum(axis=1)/(msa_len-1) #DCA_array_full.sum(axis=1)
    DCA_array_full_colSum=DCA_array_full.sum(axis=0)/(msa_len-1) #DCA_array_full.sum(axis=0)



    DCA_array_APCVar=np.zeros((DCA_array_full.shape[0],DCA_array_full.shape[1]))

    for i in range(msa_len):
        for j in range(i+1,msa_len):#
            DCA_array_APCVar[i,j]=DCA_array_full[i,j]-(DCA_array_full_rowSum[i])*(DCA_array_full_colSum[j])/(DCA_array_full_sum) 
    

    
    np.savez_compressed(DCA_folder+QueryPro_pro1+"and"+QueryPro_pro2+"_pydcaFN_array.npz",np.round(DCA_array,4))
    np.savez_compressed(DCA_folder+QueryPro_pro1+"and"+QueryPro_pro2+"_pydcaFNAPC_array.npz",np.round(DCA_array_APCVar,4))
    #return(DCA_array_APCVar) # returna uper triangular matrix 

    
# this function need conda py37_pydca enviroment 
def pydca_mfdca_FN_compresse(pp):
    QueryPro_pro1,QueryPro_pro2,L1,L2,final_MSA_folder,DCA_folder=pp
    L1,L2=int(L1),int(L2)
    
#     # this check can be avoid now as we check this before this functin to avoid always trasmite large arument and check it many times 
#     if DCA_folder+QueryPro_pro1+"and"+QueryPro_pro2+"_pydcaFNAPC_array.npz" in existed_pydcaFNAPC_files:
#         print("this output already existed")
        
#     else:
    #msa_file=pdb_benchmark_msa_folder+QueryPro_pro1+"_"+QueryPro_pro2+".fas"
    msa_file=final_MSA_folder+QueryPro_pro1+"and"+QueryPro_pro2+".fasta"


    mfdca_inst = meanfield_dca.MeanFieldDCA(
    msa_file,
    'protein',
    pseudocount = 0.5,
    seqid = 0.8,
    )

    mfdca_FN = mfdca_inst.compute_sorted_FN()
    # here in latest pydca, return value is not a matrix 
    post_pydcaFN2array_APC(QueryPro_pro1,QueryPro_pro2,L1,L2,np.asarray(mfdca_FN,dtype=object),DCA_folder)

    #time.sleep(1)# add this to prevent broken pipel problme ??
    
    return (None)



def tripleMSA_post_pydcaFN2array_APC(QueryPro_pro1,QueryPro_pro2,QueryPro_pro3,L1,L2,L3,DCA_ori,DCA_folder):

      
    DCA_array=np.zeros((L1+L2+L3,L1+L2+L3))

    array_rowIdx=np.array([t[0] for t in DCA_ori[:,0]])
    array_colIdx=np.array([t[1] for t in DCA_ori[:,0]])
    array_values=DCA_ori[:,1]
    DCA_array[array_rowIdx,array_colIdx]=array_values

    msa_len=L1+L2+L3
    DCA_array_full = DCA_array.T + DCA_array
    np.fill_diagonal(DCA_array_full,DCA_array.diagonal())# as MI of column itself is not  0!!!!


    DCA_array_full_sum=DCA_array_full.sum()/ (msa_len*msa_len-msa_len)   #DCA_array_full.sum()
    DCA_array_full_rowSum=DCA_array_full.sum(axis=1)/(msa_len-1) #DCA_array_full.sum(axis=1)
    DCA_array_full_colSum=DCA_array_full.sum(axis=0)/(msa_len-1) #DCA_array_full.sum(axis=0)


    DCA_array_APCVar=np.zeros((DCA_array_full.shape[0],DCA_array_full.shape[1]))

    for i in range(msa_len):
        for j in range(i+1,msa_len):#
            DCA_array_APCVar[i,j]=DCA_array_full[i,j]-(DCA_array_full_rowSum[i])*(DCA_array_full_colSum[j])/(DCA_array_full_sum) 
    

    
    np.savez_compressed(DCA_folder+QueryPro_pro1+"and"+QueryPro_pro2+"and"+QueryPro_pro3+"_pydcaFN_array.npz",np.round(DCA_array,4))
    np.savez_compressed(DCA_folder+QueryPro_pro1+"and"+QueryPro_pro2+"and"+QueryPro_pro3+"_pydcaFNAPC_array.npz",np.round(DCA_array_APCVar,4))
    #return(DCA_array_APCVar) # returna uper triangular matrix 

    
# this function need conda py37_pydca enviroment 
def tripleMSA_pydca_mfdca_FN_compresse(pp):
    QueryPro_pro1,QueryPro_pro2,QueryPro_pro3,L1,L2,L3,final_MSA_folder,DCA_folder=pp
    
    msa_file=final_MSA_folder+QueryPro_pro1+"and"+QueryPro_pro2+"and"+QueryPro_pro3+".fasta"
    print(msa_file)
    mfdca_inst = meanfield_dca.MeanFieldDCA(
    msa_file,
    'protein',
    pseudocount = 0.5,
    seqid = 0.8,
    )

    mfdca_FN = mfdca_inst.compute_sorted_FN()
    
    tripleMSA_post_pydcaFN2array_APC(QueryPro_pro1,QueryPro_pro2,QueryPro_pro3,L1,L2,L3,np.asarray(mfdca_FN,dtype=object),DCA_folder)

    
    return (None)








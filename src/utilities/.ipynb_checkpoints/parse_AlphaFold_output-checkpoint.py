import numpy as np
import os 
import scipy
import pickle

import json

from processANDvisulise_PDB import read_allres_pdbformat_data
from processANDvisulise_PDB import generate_proximity_matrix



def get_interprotein_cm_from_pdb(l): 
    prefix,suffix,saved_suffix,p1,p2,len1,len2=l
    pdb_file=prefix+p1+"and"+p2+suffix
    
    processedCMnpz_file=prefix+p1+"and"+p2+saved_suffix
    if not os.path.exists(processedCMnpz_file):
        pdb_struc_allRes=read_allres_pdbformat_data(pdb_file)
        contact_map,adjacency_matrix=generate_proximity_matrix(seq_1=pdb_struc_allRes,
                                                           seq_2=pdb_struc_allRes,
                                                          angstroms=10,
                                                          )
        np.savez_compressed(processedCMnpz_file,
                             contact_map) 
    else:
        try:
            contact_map=np.load(processedCMnpz_file)['arr_0'] # dont no why but sometime  np.savez_compressed in s3it mutliple processing corrupt data 
        except:
            pdb_struc_allRes=read_allres_pdbformat_data(pdb_file)
            contact_map,adjacency_matrix=generate_proximity_matrix(seq_1=pdb_struc_allRes,
                                                               seq_2=pdb_struc_allRes,
                                                              angstroms=10,
                                                              )
            np.savez_compressed(processedCMnpz_file,
                                 contact_map) 
            
    interprotein_contact_map=contact_map[0:len1,len1:]
    return(p1,p2,interprotein_contact_map)

def get_interprotein_mincm_from_pdb(l):
    p1,p2,interprotein_contact_map=get_interprotein_cm_from_pdb(l)
    return(p1,p2,interprotein_contact_map.min() )

def get_interprotein_contactprob_from_alphafold_distogram(l): 
    suffix, prefix ,p1,p2,len1,len2=l
    prediction_result_dict_file=suffix+p1+"and"+p2+prefix
    
    with open(prediction_result_dict_file, 'rb') as handle:
        prediction_result_dict = pickle.load(handle)

    prediction_result=prediction_result_dict['model_3']

    pdist = prediction_result['distogram']['logits']
    pdist = scipy.special.softmax(pdist, axis=-1)
    prob12 = np.sum(pdist[:,:,:32], axis=-1)
    interprotein_prob12=prob12[0:len1,len1:]
    
    return(p1,p2,interprotein_prob12)



def get_interprotein_contactprob_from_processedalphafoldnpz(l): 
    # input data get from test_alphafold/scipt_ParalleTest_EggNogGenereatedCustomizedMSA_multiplecomplex_deimos.py
    
    prefix,suffix,p1,p2,len1,len2=l
    contact12_npz_file=prefix+p1+"and"+p2+suffix
    try:
        prob12=np.load(contact12_npz_file)['arr_0']
        interprotein_prob12=prob12[0:len1,len1:]
    except:
        interprotein_prob12=None
    return(p1,p2,interprotein_prob12)

def get_interprotein_maxcontactprob_from_processedalphafoldnpz(l): 
    p1,p2,interprotein_prob12=get_interprotein_contactprob_from_processedalphafoldnpz(l)
    if interprotein_prob12 is not None:
        return(p1,p2,interprotein_prob12.max())
    else:
        return(p1,p2,interprotein_prob12)

def simple_precision(data,label,threshold):
    predicted_count=sum(data[label]>threshold)
    TP=sum((data[label]>threshold)&(data["status"]=="P"))
    print(f"TP:{TP};predict_P:{predicted_count};precision:{TP/predicted_count}")
    
def get_ptm_from_scoreFile(l): 
    colab_outputPath,p1,p2=l
    score_file = f"{colab_outputPath}{p1}and{p2}.custom_unrelaxed_rank_1_model_3_scores.json"
    with open(score_file) as json_file:
        pdb_scores = json.load(json_file)
    ptm=pdb_scores["ptm"]

    return(p1,p2,ptm)
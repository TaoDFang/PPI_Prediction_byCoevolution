import pandas as pd
import copy 
import random

def get_STRING1105_pdb_interact(pdb_interact_filename="/mnt/mnemo6/damian/STRING_derived_v11.5/pdb/pdb_interact.tsv",
                                retrieve_spe=511145,
                               remove_PDBoverlap=False,
                               return_overlap=False):
    
    def get_pps_withcertaintype(pdb_interact_full_spe,typelabel="complex_contact"):
        
        pdb_type=pdb_interact_full_spe.loc[pdb_interact_full_spe.iloc[:,0]==typelabel,:]
        pdb_type_PPs=pdb_type.loc[:,[2,3]].values.tolist()
        pdb_type_PPs=[tuple(sorted(tuple(pp))) for pp in pdb_type_PPs ]
        pdb_type_PPs=[("511145."+p1,"511145."+p2) for p1, p2 in pdb_type_PPs]
        pdb_type_PPs=sorted(set(pdb_type_PPs))
        
        return(pdb_type_PPs)

    pdb_interact_full=pd.read_csv(pdb_interact_filename,
                             sep="\t",header=None,index_col=None)
    pdb_interact_full_spe=pdb_interact_full.loc[pdb_interact_full.iloc[:,1]==retrieve_spe,:]


    
    pdb_interact_PPs=get_pps_withcertaintype(pdb_interact_full_spe,typelabel="complex_contact")
    pdb_complex_PPs=get_pps_withcertaintype(pdb_interact_full_spe,typelabel="complex")
    
    # here once problme is now there is overlap between  pdb_complex_PPs and pdb_interact_PPs
    # clearn pdb_complex_PPs
    print("len of intersection:",len(set(pdb_complex_PPs).intersection(set(pdb_interact_PPs))))
    
    pdb_complex_removeOverlap_PPs=list(set(pdb_complex_PPs)-set(pdb_interact_PPs))
    pdb_complexAndinteract_PPs=list(set(pdb_complex_PPs).intersection(set(pdb_interact_PPs)))
    
    if remove_PDBoverlap:
        return(list(set(pdb_interact_PPs)-set(pdb_complex_PPs)))
    else:
        if return_overlap:
            return(pdb_interact_PPs,pdb_complex_removeOverlap_PPs,pdb_complexAndinteract_PPs)
        else:
            return(pdb_interact_PPs,pdb_complex_removeOverlap_PPs)


def get_filter_pdbBenchmark(input_dataframe,
                                      pdb_interact_PPs,pdb_complex_PPs,
                           pos_lable="pdb_contact",
                           return_idx=False,
                            downsample_Ratio=False
                           ):
    metapdb_frame=copy.deepcopy(input_dataframe)
    metapdb_PPs= metapdb_frame.values.tolist()

    meta_status=list()
    for l in metapdb_PPs:
        if tuple(l[0:2]) in pdb_interact_PPs:
            meta_status.append("pdb_contact")
        elif tuple(l[0:2]) in pdb_complex_PPs:
            meta_status.append("pdb_complex")
        elif l[2]=="P":
            meta_status.append("other_pos")
        else:
            meta_status.append("N")  
    metapdb_frame["pdb_status"]=meta_status
    
    pos_idx=[i for i,s in enumerate(metapdb_frame.pdb_status) if (s==pos_lable)]
    
    if downsample_Ratio:
        neg_idx=[i for i,s in enumerate(metapdb_frame.pdb_status) if(s=="N")]
        random.seed(10)
        neg_idx=random.sample(neg_idx,downsample_Ratio)
        
    else:
        neg_idx=[i for i,s in enumerate(metapdb_frame.pdb_status) if (s=="N")]
        
        
    print("len(pos_idx),len(neg_idx)",len(pos_idx),len(neg_idx))
    sel_idx=pos_idx+neg_idx
    
    
    print(metapdb_frame.shape,len(sel_idx))
    
    if return_idx:
        return(metapdb_frame.iloc[sel_idx,:],sel_idx)
    else:
        return(metapdb_frame.iloc[sel_idx,:])
    
def add_pdbstatus(input_dataframe,pdb_interact_PPs,pdb_complex_PPs,status_colIdx=2):
    metapdb_frame=copy.deepcopy(input_dataframe)
    metapdb_PPs= metapdb_frame.values.tolist()

    meta_status=list()
    for l in metapdb_PPs:
        if tuple(l[0:2]) in pdb_interact_PPs:
            # meta_status.append("pdb_contact")
            meta_status.append("pdb_direct")
        elif tuple(l[0:2]) in pdb_complex_PPs:
            # meta_status.append("pdb_complex")
            meta_status.append("pdb_mediated")
        elif l[status_colIdx]=="P":
            meta_status.append("other_pos")
        else:
            meta_status.append("N")
            
    metapdb_frame["pdb_status"]=meta_status
    
    return(metapdb_frame)

def add_pdbstatus_fromPosFrame(input_dataframe,pdb_interact_PPs,pdb_complex_PPs):
    metapdb_frame=copy.deepcopy(input_dataframe)
    metapdb_PPs= metapdb_frame.values.tolist()

    meta_status=list()
    for l in metapdb_PPs:
        if tuple(l[0:2]) in pdb_interact_PPs:
            meta_status.append("pdb_contact")
        elif tuple(l[0:2]) in pdb_complex_PPs:
            meta_status.append("pdb_complex")
        else:
            meta_status.append("other_pos")
            
    metapdb_frame["pdb_status"]=meta_status
    
    return(metapdb_frame)

def add_pdbstatus_fromPureFrame(input_dataframe,pdb_interact_PPs,pdb_complex_PPs,):
    metapdb_frame=copy.deepcopy(input_dataframe)
    metapdb_PPs= metapdb_frame.values.tolist()

    meta_status=list()
    for l in metapdb_PPs:
        if tuple(l[0:2]) in pdb_interact_PPs:
            meta_status.append("pdb_contact")
        elif tuple(l[0:2]) in pdb_complex_PPs:
            meta_status.append("pdb_complex")    
        else:
            meta_status.append("N")
            
    metapdb_frame["pdb_status"]=meta_status
    
    return(metapdb_frame)
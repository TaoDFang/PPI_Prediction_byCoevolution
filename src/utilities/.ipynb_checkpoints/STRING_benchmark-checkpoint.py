import pandas as pd 
import copy 

def get_STRING1105_physical_interact(data_folder="",
                                retrieve_spe=511145,
                                     combined_score_thres=500,
                                     inputFrameIsReversed=True,
                                     return_dict=True,
                                    ):
    """
    here input frame shoud contain reversed version of PPI, i.e. (P1,P2) and (P2,P1)
    so is return frame/dict
    """
    STRINGcurrentSpePhyPPI_benchmark_file=f"{data_folder}{str(retrieve_spe)}.protein.physical.links.v11.5.txt.gz"
    STRINGcurrentSpePhyPPI_benchmark=pd.read_csv(STRINGcurrentSpePhyPPI_benchmark_file,sep=" ",
                                   header=0,index_col=None)

    print("STRINGcurrentSpePhyPPI_benchmark.shape:",STRINGcurrentSpePhyPPI_benchmark.shape)

    STRINGcurrentSpePhyPPI_benchmark=STRINGcurrentSpePhyPPI_benchmark.loc[STRINGcurrentSpePhyPPI_benchmark['combined_score']>combined_score_thres,:]
    print("STRINGcurrentSpePhyPPI_benchmark.shape:",STRINGcurrentSpePhyPPI_benchmark.shape)

    STRINGcurrentSpePhyPPI_benchmark=STRINGcurrentSpePhyPPI_benchmark.sort_values(by="combined_score",ascending=False)

    STRINGcurrentSpePhyPPI_benchmark.head(n=3)


    # notice here  already contains reversed version ppi 
    if return_dict==True:
        currentSpe_STRINGcurrentSpePhyPPI_posPPI=STRINGcurrentSpePhyPPI_benchmark.loc[:,["protein1","protein2"]].values.tolist()
        if inputFrameIsReversed:
            currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict={(p1,p2):1 for p1, p2 in currentSpe_STRINGcurrentSpePhyPPI_posPPI}
        else:
            currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict={tuple(sorted((p1,p2))):1 for p1, p2 in currentSpe_STRINGcurrentSpePhyPPI_posPPI}
        return(currentSpe_STRINGcurrentSpePhyPPI_posPPI_dict)
    else:
        return(STRINGcurrentSpePhyPPI_benchmark)
    
    

def add_STRINGPhyPPI_status(input_dataframe,STRINGcurrentSpePhyPPI_posPPI_dict,):
    metapdb_frame=copy.deepcopy(input_dataframe)
    metapdb_PPs= metapdb_frame.values.tolist()

    meta_status=list()
    for l in metapdb_PPs:
        if tuple(l[0:2]) in STRINGcurrentSpePhyPPI_posPPI_dict:
            meta_status.append("P")
        else:
            meta_status.append("N")
            
    metapdb_frame["STRINGPhy_status"]=meta_status
    return(metapdb_frame)


def get_string_score_dict(file_name="511145.protein.links.detailed.v11.5.txt.gz",
                         column_name_list=None):

    string_score=pd.read_csv(file_name,
                                   header=0,index_col=None,sep=" ")
    if column_name_list is None:
        string_column_names=string_score.columns.values.tolist()[2:]
    else:
        string_column_names=column_name_list
    
    string_score_dict=dict()

    for name in string_column_names:
        string_score_list=string_score.loc[:,["protein1","protein2",name]].values.tolist()
        string_score_dict[name]=dict([((p1,p2),s)for p1, p2, s in string_score_list])

    return(string_score_dict)
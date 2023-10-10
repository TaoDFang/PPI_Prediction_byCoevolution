import matplotlib.pyplot as plt 

from matplotlib.colors import LogNorm
import numpy as np 

import networkx as nx

# *********** density plot  ***********
# https://www.data-to-viz.com/graph/density2d.html
# or scatter plot with adjust marker size and transparency
def plot_list2list_scatter(msaStats_dict_ofDict,allPPs,input_tuple1,input_tuple2,feature_id,
                          plot_type="scatter",return_medianValue=True,legend=None):
    list1=[msaStats_dict_ofDict[input_tuple1][pp][feature_id] for pp in allPPs]
    list2=[msaStats_dict_ofDict[input_tuple2][pp][feature_id] for pp in allPPs]

    if return_medianValue:
        print("medianValue",np.median(list1),np.median(list2))
    if plot_type == "scatter":
        plt.scatter(list1,list2,s=1,alpha=0.8,label=legend)
        
    elif plot_type == "hist2d":
        nbins = 50
        
        fig, ax = plt.subplots()
        h=plt.hist2d(list1,list2, bins=nbins, 
                  cmap=plt.cm.BuGn_r,
                   norm=LogNorm(),
                  )
        fig.colorbar(h[3], ax=ax)
        
    plt.axline((0, 0), slope=1, color='r')
    plt.xlabel(str(input_tuple1))
    plt.ylabel(str(input_tuple2))
    #plt.show()

def plot_list2list_twofeatures_scatter(msaStats_dict_ofDict,allPPs,input_tuple,feature1_id,feature2_id):
    fea1=[msaStats_dict_ofDict[input_tuple][pp][feature1_id] for pp in allPPs]
    fea2=[msaStats_dict_ofDict[input_tuple][pp][feature2_id] for pp in allPPs]
    print(len(fea1),max(fea1))
    plt.scatter(fea1,fea2)
    plt.axline((0, 0), slope=1, color='r',alpha=0.5)
    plt.title(input_tuple)
    plt.show()
    
    

def check_PPIDist_inGraph(PPI_list,graph_PPList,AF_maxprob12_dict=None,
                          color1=None,color2=None,
                          upSel=True,maxprob12_thres=None,showRatio=True,
                         ):
    print("pp list len before filter by maxprob12",len(PPI_list))
    if AF_maxprob12_dict is not None:
        if upSel:
            PPI_list=[pp for pp in PPI_list if AF_maxprob12_dict[pp]>=maxprob12_thres]
        else:
            PPI_list=[pp for pp in PPI_list if AF_maxprob12_dict[pp]<=maxprob12_thres]
    print("pp list len after filter by maxprob12",len(PPI_list))        
        
            
    graph_PPDict={(p1,p2):1 for p1,p2  in graph_PPList}
    nx_G=nx.from_edgelist([(p1,p2) for p1,p2  in graph_PPList])


    spl = dict(nx.all_pairs_shortest_path_length(nx_G))
    dist_extend_list=list()
    for node, dis_dict in spl.items():
        dist_list=dis_dict.values()
        dist_list=[d for d in dist_list if d>0]
        dist_extend_list.extend(dist_list)
    print("mean graph dist",len(dist_extend_list),np.mean(dist_extend_list))


    coreComplex_novel_count=0
    coreComplex_nonNovel_count=0
    coreComplex_distDict=dict()
    for p1, p2 in PPI_list:

        if (p1 in spl) and (p2 in spl[p1]):
            if (p1,p2) in  graph_PPDict:
                coreComplex_nonNovel_count +=1
            else:
                coreComplex_novel_count +=1
            coreComplex_distDict[(p1,p2)]=spl[p1][p2]

    print(len(PPI_list))
    print("nonNovel,novel",coreComplex_nonNovel_count,coreComplex_novel_count)


    print(sum([1 for k, v in coreComplex_distDict.items() if v==1]))
    print(sum([1 for k, v in coreComplex_distDict.items() if v>1]))

    coreComplex_distList=[v for k, v in coreComplex_distDict.items()]
    print("mean coreComplex dist",len(coreComplex_distList),np.mean(coreComplex_distList))
    
    if showRatio:
        
        plt.hist(dist_extend_list,weights=np.ones(len(dist_extend_list)) / len(dist_extend_list),
                 bins=30,label="Full-network PP dist",
                color=color1)
        plt.hist(coreComplex_distList,weights=np.ones(len(coreComplex_distList)) / len(coreComplex_distList),
                 bins=30,label="Predicted PP dist",
                color=color2)

    else:
        plt.hist(dist_extend_list,bins=30,label="Full-network PP dist",
                color=color1)
        plt.hist(coreComplex_distList,bins=30,label="Predicted PP dist",
                color=color2)
        plt.yscale("log")
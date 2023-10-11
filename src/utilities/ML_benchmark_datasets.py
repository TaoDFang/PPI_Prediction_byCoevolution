import matplotlib.pyplot as plt
from collections import defaultdict

from PPI_benchmark_filter import get_PPIwithLimitedProFreByOr

# full STRINGPhyBalancePhyla dataset 
def return_STRINGFulllPhyBalancePhylaPPI_benchmark_PPs(ML_pos_benchmarkFrame_dict,
                                                       input_allPPs_list,
                                                       Ecoli_string_score_dict,
                                                       limitedProFre=100000,
                                                      ):
    
    currentSpe_STRINGFulllPhyBalancePhylaPPI_posPPI=[pp for pp in input_allPPs_list if pp in ML_pos_benchmarkFrame_dict]

    print("len(currentSpe_STRINGFulllPhyBalancePhylaPPI_posPPI):",len(currentSpe_STRINGFulllPhyBalancePhylaPPI_posPPI))


    posPPI_allProtins_dict={p:1 for pp in currentSpe_STRINGFulllPhyBalancePhylaPPI_posPPI for p in pp}
    print(len(posPPI_allProtins_dict))


    currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI=[pp for pp in input_allPPs_list if (pp not in ML_pos_benchmarkFrame_dict) ]

    print("len(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI):",len(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI))



    # here filter  negative proteins by only keep these pps whose bothe proteins are also in pos benchmark 

    currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI=[(p1,p2) for p1,p2  in currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI if (p1 in posPPI_allProtins_dict) and (p2 in posPPI_allProtins_dict) ]

    print("len(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI):",len(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI))
    
    # #this step can be omited if we plan to choose all neg pps 
    # random.Random(10).shuffle(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI)
    # print(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI[0:3])
    
    # constructure final positive and negigeve PPI using string score 
    currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI=[l for l in currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI if  ((l[0],l[1]) not in Ecoli_string_score_dict) ]
    print("string filtering,len(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI)",len(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI))

    ## here we  limited frquence of protein in postive and negative dataset. 
    ## when want to use all ppi, use a very large frencey number here 
    ## acturally only for negative now 
    ## futer post-filter step can be applied later , for example by Remove Ribosome protein or not 
    currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI=get_PPIwithLimitedProFreByOr(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI,limitedProFre=limitedProFre)
    print("frequency filtering ,len(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI):",len(currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI))




    # show protein frequceny distribution 
    currentSpe_STRINGPhyBalancePhylaPPI_negPPI_ProFreDict=defaultdict(int)
    for p1,p2 in currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI:
        currentSpe_STRINGPhyBalancePhylaPPI_negPPI_ProFreDict[p1] +=1
        currentSpe_STRINGPhyBalancePhylaPPI_negPPI_ProFreDict[p2] +=1
    print("len(currentSpe_STRINGPhyBalancePhylaPPI_negPPI_ProFreDict):",len(currentSpe_STRINGPhyBalancePhylaPPI_negPPI_ProFreDict))
    len([v for v in currentSpe_STRINGPhyBalancePhylaPPI_negPPI_ProFreDict.values() if v >50])
    plt.hist(currentSpe_STRINGPhyBalancePhylaPPI_negPPI_ProFreDict.values(),bins=100)
    #plt.xscale('log')
    plt.yscale('log')
    plt.show()


    return(currentSpe_STRINGFulllPhyBalancePhylaPPI_posPPI,currentSpe_STRINGFulllPhyBalancePhylaPPI_negPPI)




    



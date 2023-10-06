# keep loading Libraries here as function here may need to use in different enviroment, which only have specific Libraries

import pandas as pd


def getMetaFrame(EggNOG_maxLevel,currentSpe_TaxID,STRING_Version,DCA_thres=0,
                             given_benchmark_folder=None,
                            benchmark_suffix="STRINPhyPPI_Benchmark/",
                             splitPosandNeg=True,sort_frame=True,
                            CoEvo_data_folder="/mnt/mnemo6/tao/PPI_Coevolution/CoEvo_data_STRING11.5/"):
    '''
    
    '''
    if given_benchmark_folder is None :
        input_root_folder=CoEvo_data_folder+currentSpe_TaxID+"_EggNOGmaxLevel"+EggNOG_maxLevel+"_eggNOGfilteredData/"
        Benchmark_folder=input_root_folder+benchmark_suffix
    else:
        Benchmark_folder=given_benchmark_folder
    
    print("Benchmark_folder:",Benchmark_folder)
    allPPI_allInfo_frame=pd.read_csv(Benchmark_folder+"allPPI_allInfo_frame.csv",
                                 header=0,index_col=None,sep="\t")
    
    print("allPPI_allInfo_frame.shape:",allPPI_allInfo_frame.shape)
    
    if splitPosandNeg:
        Predictive_DCA=DCA_thres

        Pos_allPPI_allInfo_frame=allPPI_allInfo_frame.loc[allPPI_allInfo_frame['benchmark_status']=="P",:]
        #print("Pos_allPPI_allInfo_frame.shape:",Pos_allPPI_allInfo_frame.shape)
        Pos_allPPI_allInfo_frame=Pos_allPPI_allInfo_frame.loc[Pos_allPPI_allInfo_frame['maxBetDCA_score']>=Predictive_DCA,:]
        print("Pos_allPPI_allInfo_frame.shape:",Pos_allPPI_allInfo_frame.shape)


        Neg_allPPI_allInfo_frame=allPPI_allInfo_frame.loc[allPPI_allInfo_frame['benchmark_status']=="N",:]
        #print("Neg_allPPI_allInfo_frame.shape:",Neg_allPPI_allInfo_frame.shape)
        Neg_allPPI_allInfo_frame=Neg_allPPI_allInfo_frame.loc[Neg_allPPI_allInfo_frame['maxBetDCA_score']>=Predictive_DCA,:]
        print("Neg_allPPI_allInfo_frame.shape:",Neg_allPPI_allInfo_frame.shape)

        if sort_frame:
            Pos_allPPI_allInfo_frame=Pos_allPPI_allInfo_frame.sort_values(by="maxBetDCA_score",ascending=False)
            Neg_allPPI_allInfo_frame=Neg_allPPI_allInfo_frame.sort_values(by="maxBetDCA_score",ascending=False)

        allPPI_allInfo_frame= pd.concat([Pos_allPPI_allInfo_frame,Neg_allPPI_allInfo_frame])
        print("allPPI_allInfo_frame.shape:",allPPI_allInfo_frame.shape)

    return(allPPI_allInfo_frame)



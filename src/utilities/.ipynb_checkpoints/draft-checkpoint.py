def get_BestHomologousDCAs_top5DCAs_fromMultiSpes(Query_allPPI_allInfo_dict,
                                    Query_allPPI_top5DCAs_dict,
                                    BestHomologousPP_Subject_allPPI_top5DCAs_listDict,
                                    Query2Subject_BestHomologous_ignoreSubjectDCA_dict_listDict,
                                        topDCA_num=5,
                                    with_status=True,
                                                 ):
    
    Query_BestHomologousDCAs_dict=defaultdict(list)

    for sub_PP in Query_allPPI_allInfo_dict:
        #print(Query_allPPI_allInfo_dict[sub_PP])
        if with_status:
            Query_BestHomologousDCAs_dict[sub_PP]=list(Query_allPPI_allInfo_dict[sub_PP][0])
        else:
            Query_BestHomologousDCAs_dict[sub_PP]=[]
            
        Query_BestHomologousDCAs_dict[sub_PP].extend(Query_allPPI_top5DCAs_dict[sub_PP][0:topDCA_num])

        #for Subject_phylum, Subject_speID in Subject_tupleList:
        for Subject_phylum, Subject_speID in Query2Subject_BestHomologous_ignoreSubjectDCA_dict_listDict.keys():
            #print("Subject_phylum, Subject_speID:",Subject_phylum, Subject_speID)
            BestHomologousPP_Subject_allPPI_top5DCAs_dict=BestHomologousPP_Subject_allPPI_top5DCAs_listDict[(Subject_phylum, Subject_speID)]
            Query2Subject_BestHomologous_ignoreSubjectDCA_dict=Query2Subject_BestHomologous_ignoreSubjectDCA_dict_listDict[(Subject_phylum, Subject_speID)]

            best_Subject_PP=Query2Subject_BestHomologous_ignoreSubjectDCA_dict[sub_PP] if sub_PP in Query2Subject_BestHomologous_ignoreSubjectDCA_dict else []
            #print(best_Subject_PP)
            if len(best_Subject_PP)==0:  # here Subject_PP is tuple with size 2
                Query_BestHomologousDCAs_dict[sub_PP].extend([np.nan for i in range(topDCA_num)])
            else:

                if best_Subject_PP in BestHomologousPP_Subject_allPPI_top5DCAs_dict:
                    Query_BestHomologousDCAs_dict[sub_PP].extend(BestHomologousPP_Subject_allPPI_top5DCAs_dict[best_Subject_PP][0:topDCA_num])

                else:
                    Query_BestHomologousDCAs_dict[sub_PP].extend([np.nan for i in range(topDCA_num)])

                        
    return(Query_BestHomologousDCAs_dict)

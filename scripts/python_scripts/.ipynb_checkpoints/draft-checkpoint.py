


def allQuery2SubjectSingleProteinMapping(Query2Subject_QueSpeAllPPI_homologous_dict_listDict,
                                                                  ):
 

    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict=defaultdict(dict)
    for phylum_speID,temp_Query2Subject_QueSpeAllPPI_homologous_dict in Query2Subject_QueSpeAllPPI_homologous_dict_listDict.items():
        print(phylum_speID)
        for Query_pp, Subject_pps in temp_Query2Subject_QueSpeAllPPI_homologous_dict.items():
            Query_pp1,Query_pp2=Query_pp
            for Subject_pp in Subject_pps:
    #             if(len(Subject_pps)==1):
    #                 print(Subject_pp)
                Subject_pp1,Subject_pp2=Subject_pp

                # here we should notice one  single protein can have multiple homologous protein pairs
                # and same one singel protein to single protein mapping can appear many times in different homologous protein pair mapping 
                if Query_pp1 not in Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID]:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID][Query_pp1]={Subject_pp1} #set
                else:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID][Query_pp1].add(Subject_pp1)

                if Query_pp2 not in Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID]:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID][Query_pp2]={Subject_pp2} # set 
                else:
                    Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict[phylum_speID][Query_pp2].add(Subject_pp2)


    return(Query2Subject_QueSpeAllPPI_homologous_singleProteinMaping_listDict)
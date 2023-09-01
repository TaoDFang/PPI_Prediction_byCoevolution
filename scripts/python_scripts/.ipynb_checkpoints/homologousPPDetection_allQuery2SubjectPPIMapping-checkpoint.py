import argparse
import os 
import glob
import pickle
import multiprocessing as mp

from IntergrateBestHomologousPPCoEvo_unNameSorted import allQuery2SubjectPPIMapping_getNameUnSorted_COGasLink

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()


    parser.add_argument('-q','--Query_tuple', type=str, help='Query_tuple')
    parser.add_argument('-s','--Subject_tupleList', type=str, help='Subject_tupleList')
    parser.add_argument('-p','--PPIInfoBeforeCoEvoComp_csv', type=str, help='PPIInfoBeforeCoEvoComp_csv')
    parser.add_argument('-t','--homologous_COG2PP_path', type=str, help='homologous_COG2PP_path')
    parser.add_argument('-m','--homologous_allQuery2SubjectPPIMapping_path', type=str, help='homologous_allQuery2SubjectPPIMapping_path')

    args = parser.parse_args()
    Query_tuple=args.Query_tuple
    Subject_tupleList=args.Subject_tupleList
    PPIInfoBeforeCoEvoComp_csv=args.PPIInfoBeforeCoEvoComp_csv
    homologous_COG2PP_path=args.homologous_COG2PP_path
    homologous_allQuery2SubjectPPIMapping_path=args.homologous_allQuery2SubjectPPIMapping_path

    Query_tuple=Query_tuple.split("_")
    Subject_tupleList=Subject_tupleList.split("_")
    Subject_tupleList=[(Subject_tupleList[i],Subject_tupleList[i+1]) for i in range(0,len(Subject_tupleList)-1,2)]
    
    print("Query_tuple:",Query_tuple)
    print("Subject_tupleList:",Subject_tupleList)
    Query2Subject_QueSpeAllPPI_homologous_ignoreQueryDCA_dict_listDict=allQuery2SubjectPPIMapping_getNameUnSorted_COGasLink(Query_tuple,Subject_tupleList,PPIInfoBeforeCoEvoComp_csv,homologous_COG2PP_path)


    with open(homologous_allQuery2SubjectPPIMapping_path+"NameUnsorted_Query2Subject_QueSpeAllPPI_homologous_dict_listDict.pickle", 'wb') as handle:
            pickle.dump(Query2Subject_QueSpeAllPPI_homologous_ignoreQueryDCA_dict_listDict,handle)

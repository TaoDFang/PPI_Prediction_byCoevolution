#!/usr/bin/env nextflow

nextflow.enable.dsl=2



params.newSTRING_rootFolder="${params.PPI_Coevolution}/STRING_data_11.5" //this folder is just newSTRING_rootFolder
params.homologous_ppPath="${params.newSTRING_rootFolder}/homologous_pp"
params.homologous_SeqMappingPath="${params.homologous_ppPath}/SeqMapping" 

params.CoEvo_data_folder="${params.PPI_Coevolution}/CoEvo_data_STRING11.5/"


include {RawFastaFilesAndMetaData_workflow } from './RawFastaFilesAndMetaData_workflow.nf'
include {prepareSingleMSA_workflow} from './prepareSingleMSA_workflow.nf'
include {prepareSingleMSA_workflow as subject1_prepareSingleMSA_workflow} from './prepareSingleMSA_workflow.nf'
include {prepareSingleMSA_workflow as subject2_prepareSingleMSA_workflow} from './prepareSingleMSA_workflow.nf'
include {prepareSingleMSA_workflow as subject3_prepareSingleMSA_workflow} from './prepareSingleMSA_workflow.nf'
include {Query_PairedMSA_preprocessing_workflow} from  "./Query_PairedMSA_preprocessing_workflow.nf"


spe_list_ch= Channel.value( ["${params.query_currentSpe_TaxID}", "${params.subject1_currentSpe_TaxID}", 
                             "${params.subject2_currentSpe_TaxID}", "${params.subject3_currentSpe_TaxID}"] )

Query_tuple_ch=Channel.value(["${params.query_current_EggNOG_maxLevel}","${params.query_currentSpe_TaxID}"])
Subject_tupleList_ch=Channel.value(["${params.subject1_current_EggNOG_maxLevel}","${params.subject1_currentSpe_TaxID}",
                                "${params.subject2_current_EggNOG_maxLevel}","${params.subject2_currentSpe_TaxID}",
                                "${params.subject3_current_EggNOG_maxLevel}","${params.subject3_currentSpe_TaxID}"])


// the main logic is from http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_homologousPPDetection.ipynb
// and/or mainly from  http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_homologousPPDetection_STRINGPhyBalancePhyla.ipynb
// and http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/Compute_allPPI_homologousPPDetection.ipynb

workflow homologousPPDetectionAndCompuation_workflow{    
    take:
        spe_list_ch
        eggNOG_folder
        Query_tuple_ch
        Subject_tupleList_ch
        query_PPIInfoBeforeCoEvoComp_csv
        SubjectProSeqPath_ByProteins_ch
        QueryProSeqPath_ByProtein_ch
        homologous_SeqMappingPath_ch
        newSTRING_rootFolder_ch
        CoEvo_data_folder_ch
        
    main:
        println "scriptFile: " + workflow.scriptFile
        println "projectDir: " + workflow.projectDir
        println "launchDir: " + workflow.launchDir
        println "workDir: " + workflow.workDir
        println "configFiles: " + workflow.configFiles
    
    
        homologousPPDetection_COG2PPMapping_ch=homologousPPDetection_COG2PPMapping(spe_list_ch,RawFastaFilesAndMetaData_workflow.out.eggNOG_folder)

        homologousPPDetection_allQuery2SubjectPPIMapping_ch=homologousPPDetection_allQuery2SubjectPPIMapping(Query_tuple_ch,Subject_tupleList_ch,
                                                                      Query_PairedMSA_preprocessing_workflow_ch.PPIInfoBeforeCoEvoComp_csv,
                                                                      homologousPPDetection_COG2PPMapping_ch.homologous_COG2PP_path)

        // here need to change to exclude query species
        // try with regular expression failed, use nextflow filter operator ? https://stackoverflow.com/questions/76104581/nextflow-regex-on-path
        SubjectProSeqPath_ByProteins_ch=Channel
            .fromPath("${params.newSTRING_rootFolder}/*ByProteins/",type:"dir")
            .filter( ~/^((?!${params.query_currentSpe_TaxID}).)*/ ) // this works ! https://stackoverflow.com/questions/406230/regular-expression-to-match-a-line-that-doesnt-contain-a-word
            // .filter( ~/.*511145.*/ )  // this works 
            // .filter( ~/.*(?!511145).*/ )  // this does not work, why ???

        // SubjectSpe_MiddleDatas_ch=Channel
        //     .fromPath("${params.newSTRING_rootFolder}/*MiddleData/",type:"dir") // here use this path to extract both phylum and species id, can alos use maps in nextflow 
        //     .filter( ~/^((?!511145).)*/ ) 


        QueryProSeqPath_ByProtein_ch=Channel
            .fromPath("${params.newSTRING_rootFolder}/*ByProteins/",type:"dir")
            .filter( ~/.*${params.query_currentSpe_TaxID}.*/ )  // this works 




        homologousPPDetection_SeqMapping_ch=homologousPPDetection_SeqMapping(Query_tuple_ch,
                                                                             Subject_tupleList_ch,
                                                                             QueryProSeqPath_ByProtein_ch,
                                                                             SubjectProSeqPath_ByProteins_ch,
                                                                            homologousPPDetection_allQuery2SubjectPPIMapping_ch.homologous_allQuery2SubjectPPIMapping_path,
                                                                            )


        homologous_SeqMappingPath_ch=Channel.fromPath("${params.homologous_SeqMappingPath}",type:"dir")
        homologousPPDetection_allQuery2SubjectPPIMapping_BestHomologous_ch=homologousPPDetection_allQuery2SubjectPPIMapping_BestHomologous(Query_tuple_ch,
                                                                                    Subject_tupleList_ch, homologous_SeqMappingPath_ch,
                                                                                    homologousPPDetection_allQuery2SubjectPPIMapping_ch.homologous_allQuery2SubjectPPIMapping_path)

        newSTRING_rootFolder_ch=Channel.fromPath("${params.newSTRING_rootFolder}",type:"dir") 
        homologousPPDetection_preparePairedMSA_ch=homologousPPDetection_preparePairedMSA(Query_tuple_ch,
                                                                                    Subject_tupleList_ch,
                                                                                    newSTRING_rootFolder_ch,
                    homologousPPDetection_allQuery2SubjectPPIMapping_BestHomologous_ch.homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path)

        CoEvo_data_folder_ch=Channel.fromPath("${params.CoEvo_data_folder}",type:"dir") 
        homologousPPDetection_ComputeHomologousDCA_ch=homologousPPDetection_ComputeHomologousDCA(Query_tuple_ch,Subject_tupleList_ch,
                                                                                        newSTRING_rootFolder_ch,
                                                                                        CoEvo_data_folder_ch)


    
}

workflow {
    
    
    // ************Download all metadata and sequencing data**********
    RawFastaFilesAndMetaData_workflow() 
    
   // ************Prepare single MSA **********
    //http://localhost:8206/lab/workspaces/auto-I/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
    //http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.py
    query_prepareSingleMSA_workflow=prepareSingleMSA_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                         Channel.value("${params.query_current_EggNOG_maxLevel}"),
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                              RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                              RawFastaFilesAndMetaData_workflow.out.species_file,
                              RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                             )
    
    // this is to prevent error  
    //Process 'prepareSingleMSA_ParseCurSpeFastaByProteins' has been already used -- If you need to reuse the same component, include it with a different name or include it in a different workflow context    
    subject1_prepareSingleMSA_workflow=subject1_prepareSingleMSA_workflow(Channel.value("${params.subject1_currentSpe_TaxID}"),
                         Channel.value("${params.subject1_current_EggNOG_maxLevel}"),
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                              RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                              RawFastaFilesAndMetaData_workflow.out.species_file,
                              RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                             )

    subject2_prepareSingleMSA_workflow=subject2_prepareSingleMSA_workflow(Channel.value("${params.subject2_currentSpe_TaxID}"),
                         Channel.value("${params.subject2_current_EggNOG_maxLevel}"),
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                              RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                              RawFastaFilesAndMetaData_workflow.out.species_file,
                              RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                             )

    subject3_prepareSingleMSA_workflow=subject3_prepareSingleMSA_workflow(Channel.value("${params.subject3_currentSpe_TaxID}"),
                         Channel.value("${params.subject3_current_EggNOG_maxLevel}"),
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                              RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                              RawFastaFilesAndMetaData_workflow.out.species_file,
                              RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                             )

    // ************Prepare paired MSA of qurey speceis, then later detect all subject homologous pp for these pp in query species**********
    // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
    // here for the because the output generate millions of output,  to prevent completely re-run of progrom in case of accidentaly interuption?
    // some input channle should be the final output path, not the the temperaoy woring directory of previous path 
    // but actully if last step finished correctly , its okay, more the problem of last step compute DCA
    Query_PairedMSA_preprocessing_workflow_ch=Query_PairedMSA_preprocessing_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                                           query_prepareSingleMSA_workflow.currentSpeMSAGapsFilteringMetaFolder,
                                           query_prepareSingleMSA_workflow.currentSpe_msa_removeGaps_path)
    

    
    
    
}




/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run query2subject_homologousPPDetection_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

*/

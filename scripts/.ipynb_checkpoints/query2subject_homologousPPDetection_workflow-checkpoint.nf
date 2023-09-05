#!/usr/bin/env nextflow

nextflow.enable.dsl=2



params.newSTRING_rootFolder="${params.PPI_Coevolution}/STRING_data_11.5" //this folder is just newSTRING_rootFolder
params.homologous_ppPath="${params.PPI_Coevolution}/STRING_data_11.5/homologous_pp" 


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
workflow {
    
    println "scriptFile: " + workflow.scriptFile
    println "projectDir: " + workflow.projectDir
    println "launchDir: " + workflow.launchDir
    println "workDir: " + workflow.workDir
    println "configFiles: " + workflow.configFiles
    
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

    // ************Prepare paired MSA **********
    // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
    // here for the because the output generate millions of output,  to prevent completely re-run of progrom in case of accidentaly interuption?
    // some input channle should be the final output path, not the the temperaoy woring directory of previous path 
    // but actully if last step finished correctly , its okay, more the problem of last step compute DCA
    Query_PairedMSA_preprocessing_workflow_ch=Query_PairedMSA_preprocessing_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                                           query_prepareSingleMSA_workflow.currentSpeMSAGapsFilteringMetaFolder,
                                           query_prepareSingleMSA_workflow.currentSpe_msa_removeGaps_path)
    

    homologousPPDetection_COG2PPMapping_ch=homologousPPDetection_COG2PPMapping(spe_list_ch,RawFastaFilesAndMetaData_workflow.out.eggNOG_folder)
    
    homologousPPDetection_COG2PPMapping_ch=homologousPPDetection_allQuery2SubjectPPIMapping(Query_tuple_ch,Subject_tupleList_ch,
                                                                  Query_PairedMSA_preprocessing_workflow_ch.PPIInfoBeforeCoEvoComp_csv,
                                                                  homologousPPDetection_COG2PPMapping_ch.homologous_COG2PP_path)
    
    // here need to change to exclude query species
    // try with regular expression failed, use nextflow filter operator ? https://stackoverflow.com/questions/76104581/nextflow-regex-on-path
    SubjectProSeqPath_ByProteins_ch=Channel
        .fromPath("${params.newSTRING_rootFolder}/*ByProteins/",type:"dir")
        .filter( ~/^((?!511145).)*/ ) // this works ! https://stackoverflow.com/questions/406230/regular-expression-to-match-a-line-that-doesnt-contain-a-word
        // .filter( ~/.*511145.*/ )  // this works 
        // .filter( ~/.*(?!511145).*/ )  // this does not work, why ???
        // .filter{ ~/.*(?!511145.)*/ }
        
    
        // .filter{it =~ /.*511145.*/}
        // .filter{ ~/.*511145.*/ }
        // .filter( ~/.*/ )
        // .filter( ".*511145.*")
        // .filter( ^(?!511145).* )
    
    // // seems regular expression works on val but not on path ??!! https://stackoverflow.com/questions/76104581/nextflow-regex-on-path
    // SubjectProSeqPath_ByProteins_ch=Channel
    //                                 .of( 'a', 'b', 'aa', 'adfafbc', 3, 4.5 )
    //                                 .filter( ~/.*fa.*/ )
    
    
    
    homologousPPDetection_SeqMapping_ch=homologousPPDetection_SeqMapping(SubjectProSeqPath_ByProteins_ch)
    
    
}







    

// prepare homoglogous COG group to pp group  mapping in name unsorted manner. so we can know single protein mapping relations
// and here we use homologou pp under eggenog level  max level 2 (root level ) 
// so information stored is indepedent on pps in query species benchmark datasete in which eggno level  could be different 
process homologousPPDetection_COG2PPMapping {
    
    publishDir "${params.homologous_ppPath}", mode: "copy"
    
    label "simple_py_process"
    
    debug true //echo true echo directive is deprecated
    
    
    input: 
        val  spe_list_ch
        path eggNOG_folder
        
    output:
        path "COG2PP/", emit: homologous_COG2PP_path
    script:
        
    """
        homologous_COG2PP_path="COG2PP/" 
        mkdir -p \${homologous_COG2PP_path}
        
        # here -s need to a array, or multiple values without "," in betwween
        # didnt figure out how to use nextflow value to bash array, so change it the intermediated string 
        #echo \${spe_list_ch}
        #echo \${spe_list_ch.toList()}
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/homologousPPDetection_COG2PPMapping.py -s ${spe_list_ch.join("_")} \
        -egg "${eggNOG_folder}/groups/" -t "\${homologous_COG2PP_path}" -n ${params.small_mp_task_nums}
    """
    
}



// do protein mapping for homologous pp and single protein mapping, not the one best homologous pp
process homologousPPDetection_allQuery2SubjectPPIMapping {
    
    publishDir "${params.newSTRING_rootFolder}", mode: "copy"
    
    
    label "simple_py_process"
    
    debug true //echo true echo directive is deprecated
    
    
    input: 
        val Query_tuple_ch
        val Subject_tupleList_ch
        path PPIInfoBeforeCoEvoComp_csv
        path homologous_COG2PP_path
        
    output:
        path "${Query_tuple_ch[1]}_EggNOGmaxLevel${Query_tuple_ch[0]}_allQuery2SubjectPPIMapping/", emit:homologous_allQuery2SubjectPPIMapping_path
    script:
        
    """
        homologous_allQuery2SubjectPPIMapping_path="${Query_tuple_ch[1]}_EggNOGmaxLevel${Query_tuple_ch[0]}_allQuery2SubjectPPIMapping/" 
        mkdir -p \${homologous_allQuery2SubjectPPIMapping_path}

        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/homologousPPDetection_allQuery2SubjectPPIMapping.py \
        -q ${Query_tuple_ch.join("_")} -s ${Subject_tupleList_ch.join("_")} \
        -p ${PPIInfoBeforeCoEvoComp_csv} -t "${homologous_COG2PP_path}/" \
        -m \${homologous_allQuery2SubjectPPIMapping_path}
        
    """
    
}


// do blastp/mapping between single query proteins and single subject proteins 
process homologousPPDetection_SeqMapping {
    
    publishDir "${params.homologous_ppPath}", mode: "copy"
    
    
    label "simple_py_process"
    
    debug true //echo true echo directive is deprecated
    
    
    input: 
        path SubjectProSeqPath_ByProtein_ch
        // val SubjectProSeqPath_ByProtein_ch
        
    output:
        
    script:
        
    """
        homologous_SeqMappingPath="SeqMapping/" 
        mkdir -p \${homologous_SeqMappingPath}
        
        echo "process homologousPPDetection_SeqMapping: ${SubjectProSeqPath_ByProtein_ch}"


        
    """
    
}

//         export PYTHONPATH="${projectDir}/../src/utilities/" 
//         python ${projectDir}/python_scripts/homologousPPDetection_SeqMapping.py





// // prepare paired MSA for homologous pp in other species 
// process homologousPPDetection_preparePairedMSA {
    
//     publishDir "${params.}", mode: "copy"
    
    
//     label "simple_py_process"
    
//     debug true //echo true echo directive is deprecated
    
    
//     input: 

//     output:
       
        
//     script:
        
//     """


//         export PYTHONPATH="${projectDir}/../src/utilities/" 
//         python ${projectDir}/python_scripts/homologousPPDetection_preparePairedMSA.py 
        
//     """
    
// }









/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run query2subject_homologousPPDetection_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

*/


// !!!!
//for need to do DCA compuation for millions of PPs , if interupt, better not re-start from scratch 
// in this case, better not use output to the temperory working directory and then copy to publishDir
// but direct use final output folder as the input channel or as input parameters 

process coevolutionComputation_mfDCA_createdFolder {
    debug true //echo true echo directive is deprecated
     // here is important to not use publishDir, otherwisethe  DCA_coevolutoin_path will get overwroten
    // output:
    script: 
    """      
    mkdir -p ${params.DCA_coevolutoin_path}
    """
}

process coevolutionComputation_mfDCA {
    // publishDir "${params.input_root_folder}",mode: "copy"
    debug true //echo true echo directive is deprecated
    
    input: 
        path DCA_coevolutoin_path
        path currentSpeMiddleDataPath
        path PPIInfoBeforeCoEvoComp_csv
        path pairedMSA_Nf90_folder
    // output:
    script: 
    """      
        
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 

        python ${projectDir}/python_scripts/coevolutionComputation_mfDCA.py -dpath "${DCA_coevolutoin_path}/"  \
        -m "${currentSpeMiddleDataPath}/" -acsv ${PPIInfoBeforeCoEvoComp_csv} -nf90f "${pairedMSA_Nf90_folder}/" \
        -n ${params.middle_mp_task_nums}
        
    """
}
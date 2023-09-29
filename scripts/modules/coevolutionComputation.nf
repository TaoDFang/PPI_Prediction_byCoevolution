
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
    echo process coevolutionComputation_mfDCA_createdFolder
    echo  $CONDA_DEFAULT_ENV
    
    mkdir -p ${params.DCA_coevolutoin_path}
    """
}

process coevolutionComputation_mfDCA {
    // publishDir "${params.input_root_folder}",mode: "copy"
    debug true //echo true echo directive is deprecated
    
    label "coevolutionComputation_mfDCA_process"
    
    input: 
        path DCA_coevolutoin_path
        path currentSpeMSAGapsFilteringMetaFolder
        path PPIInfoBeforeCoEvoComp_csv
        path pairedMSA_Nf90_folder
    // output:
    script: 
    """      
        conda info --envs 
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
    
        python ${projectDir}/python_scripts/coevolutionComputation_mfDCA.py -dpath "${DCA_coevolutoin_path}/"  \
        -m "${currentSpeMSAGapsFilteringMetaFolder}/" -acsv ${PPIInfoBeforeCoEvoComp_csv} -nf90f "${pairedMSA_Nf90_folder}/" \
        -n ${params.middle_mp_task_nums}
        
    """
}
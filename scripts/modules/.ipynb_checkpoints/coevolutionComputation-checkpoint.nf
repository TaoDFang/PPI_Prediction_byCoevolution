process coevolutionComputation_mfDCA_createdFolder {
    debug true //echo true echo directive is deprecated
    label "simple_py_process"
     // here is important to not use publishDir, otherwisethe  DCA_coevolutoin_path will get overwroten
    // output:
    script: 
    """      
    echo process coevolutionComputation_mfDCA_createdFolder
    echo  $CONDA_DEFAULT_ENV
    
    mkdir -p ${params.DCA_coevolutoin_path}
   
    
    """
}
 // mkdir -p ${params.IndexDCA_coevolutoin_path}

process coevolutionComputation_mfDCA_preparaIndexFile {
    publishDir "${params.input_root_folder}",mode: "copy"
    debug true //echo true echo directive is deprecated
    
    label "simple_py_process"
    
    input: 
        path DCA_coevolutoin_path  
        path currentSpeMSAGapsFilteringMetaFolder
        path PPIInfoBeforeCoEvoComp_csv
        path pairedMSA_Nf90_folder
        val DCA_blockNum_ch
    output:
        path "coevolutoin_computation_IndexDCA/",type: "dir", emit: IndexDCA_coevolutoin_path
    script: 
    """      
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        
        echo "DCA_coevolutoin_path:${DCA_coevolutoin_path}"  # here we want to get full path for next steps 
        echo "DCA_coevolutoin_path:${pairedMSA_Nf90_folder}"


        DCA_coevolutoin_path_full=\$(realpath "${DCA_coevolutoin_path}")
        echo "DCA_coevolutoin_path_full:\${DCA_coevolutoin_path_full}"
        
        pairedMSA_Nf90_folder_full=\$(realpath "${pairedMSA_Nf90_folder}")
        echo "pairedMSA_Nf90_folder_full:\${pairedMSA_Nf90_folder_full}"


        
        
        IndexDCA_coevolutoin_path="coevolutoin_computation_IndexDCA/"
        mkdir -p \${IndexDCA_coevolutoin_path}
        

        python ${projectDir}/python_scripts/coevolutionComputation_mfDCA_preparaIndexFile.py  \
        -dpath "\${DCA_coevolutoin_path_full}/"  -dipath "\${IndexDCA_coevolutoin_path}"    \
        -m "${currentSpeMSAGapsFilteringMetaFolder}/" -acsv ${PPIInfoBeforeCoEvoComp_csv} -nf90f "\${pairedMSA_Nf90_folder_full}/" \
        -bn  ${DCA_blockNum_ch} -n ${params.middle_mp_task_nums}
    

        
    """
}





// force this process to run ofterwards, have to !!!!!!
process coevolutionComputation_mfDCA_parallel {
    debug true //echo true echo directive is deprecated
    
    label "coevolutionComputation_mfDCA_parallel_process"
    
    input: 
        path IndexDCA_coevolutoin_path
        each idxCH
    // output:
    script: 
    """      
        conda info --envs 
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        
        echo "idxCH:${idxCH}"
    
        python ${projectDir}/python_scripts/coevolutionComputation_mfDCA_parallel.py -dipath "${IndexDCA_coevolutoin_path}/" \
        -i "${idxCH}"
    

        
    """
}






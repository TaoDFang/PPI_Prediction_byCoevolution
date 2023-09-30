
// !!!!
//for need to do DCA compuation for millions of PPs , if interupt, better not re-start from scratch 
// in this case, better not use output to the temperory working directory and then copy to publishDir
// but direct use final output folder as the input channel or as input parameters 

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
    // adapt from  http://localhost:8206/lab/tree/code/MNF/notebooks/ScienceCluster_code/STRING_Data_11.5/Compute_allPPI.ipynb
    // publishDir "${params.input_root_folder}",mode: "copy"
    debug true //echo true echo directive is deprecated
    
    label "simple_py_process"
    
    input: 
        path DCA_coevolutoin_path  // one question is if input will triper current process is the its content changed?,  Channel.fromPath(params.DCA_coevolutoin_path,type:'dir')
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
        
        #notice here outsie it ( not {, and differentation between nextflow and bash variabl
        #https://itslinuxfoss.com/fix-bash-bad-subtituion-linux/#:~:text=The%20error%20%E2%80%9Cbad%20substitution%E2%80%9D%20is%20the%20wrong%20syntax%2C%20repetition,spaces%20from%20the%20bash%20script.
        #https://github.com/nextflow-io/nextflow/discussions/2887

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








// force this process to run ofter this, have to !!!!!!
process coevolutionComputation_mfDCA_parallel {
    // adapt from  http://localhost:8206/lab/tree/code/MNF/notebooks/ScienceCluster_code/STRING_Data_11.5/Compute_allPPI.ipynb
    // publishDir "${params.input_root_folder}",mode: "copy"
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





// process coevolutionComputation_mfDCA {
//     // publishDir "${params.input_root_folder}",mode: "copy"
//     debug true //echo true echo directive is deprecated
    
//     label "coevolutionComputation_mfDCA_process"
    
//     input: 
//         path DCA_coevolutoin_path
//         path currentSpeMSAGapsFilteringMetaFolder
//         path PPIInfoBeforeCoEvoComp_csv
//         path pairedMSA_Nf90_folder
//     // output:
//     script: 
//     """      
//         conda info --envs 
        
//         export PYTHONPATH="${projectDir}/../src/utilities/" 
    
//         python ${projectDir}/python_scripts/coevolutionComputation_mfDCA.py -dpath "${DCA_coevolutoin_path}/"  \
//         -m "${currentSpeMSAGapsFilteringMetaFolder}/" -acsv ${PPIInfoBeforeCoEvoComp_csv} -nf90f "${pairedMSA_Nf90_folder}/" \
//         -n ${params.middle_mp_task_nums}
        
//     """
// }





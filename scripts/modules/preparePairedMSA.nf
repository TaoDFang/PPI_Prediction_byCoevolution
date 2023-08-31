process preparePairedMSA_oneRunWithNf90 {
    publishDir "${params.input_root_folder}",mode: "copy"
    
    label "many_cpu_process"
    
    debug true //echo true echo directive is deprecated
    
    input: 
        val currentSpe_TaxID_ch
        path currentSpeMSAGapsFilteringMetaFolder
        path currentSpe_msa_removeGaps_path

    output:
        path "pair_MSA_unfiltered_PasteAlign/", emit: pairedMSA_unfiltered_folder
        path "pair_MSA_hhfilter_PasteAlign/", emit:pairedMSA_hhfilter_folder
        path "pair_MSA_Nf90_PasteAlign/", emit:pairedMSA_Nf90_folder
        path "pair_MSA_Nf90.csv", emit:pairedMSA_Nf90_csv
        
    script: 
    """      
        pairedMSA_unfiltered_folder="pair_MSA_unfiltered_PasteAlign/"
        pairedMSA_hhfilter_folder="pair_MSA_hhfilter_PasteAlign/"
        pairedMSA_Nf90_folder="pair_MSA_Nf90_PasteAlign/"
        pairedMSA_Nf90_csv="pair_MSA_Nf90.csv"

        
        mkdir -p \${pairedMSA_unfiltered_folder}
        mkdir -p \${pairedMSA_hhfilter_folder}
        mkdir -p \${pairedMSA_Nf90_folder}

        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/preparePairedMSA_oneRunWithNf90.py -i ${currentSpe_TaxID_ch} \
        -m "${currentSpeMSAGapsFilteringMetaFolder}/" -rp "${currentSpe_msa_removeGaps_path}/" \
        -un \${pairedMSA_unfiltered_folder} -hh \${pairedMSA_hhfilter_folder} \
        -nf90f \${pairedMSA_Nf90_folder} -nf90csv \${pairedMSA_Nf90_csv} \
        -nf90 ${params.Nf90_thres} -n ${params.large_mp_task_nums}
        
    """
}


process preparePairedMSA_removeHomologousPairs {
    publishDir "${params.input_root_folder}",mode: "copy"
    
    label "many_cpu_process"
    
    debug true //echo true echo directive is deprecated
    
    input: 
        path pairedMSA_Nf90_csv
        path pairedMSA_Nf90_folder

    output:
        path "sameProteinRatio.csv", emit: pairedMSA_sameProteinRatio_csv
        path "PPIInfoBeforeCoEvoComp.csv", emit: PPIInfoBeforeCoEvoComp_csv
        
    script: 
    """      
        pairedMSA_sameProteinRatio_csv="sameProteinRatio.csv" # no space around "=" !!
        PPIInfoBeforeCoEvoComp_csv="PPIInfoBeforeCoEvoComp.csv"

        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/preparePairedMSA_removeHomologousPairs.py \
        -nf90csv ${pairedMSA_Nf90_csv} -nf90f "${pairedMSA_Nf90_folder}/" \
        -scsv \${pairedMSA_sameProteinRatio_csv} -acsv \${PPIInfoBeforeCoEvoComp_csv=} \
        -n ${params.large_mp_task_nums}  
    """
}
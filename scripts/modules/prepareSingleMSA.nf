
// ************Prepare single MSA **********

//section ParseEcoliFastaByProtein, need in phmmer section to 
process prepareSingleMSA_ParseCurSpeFastaByProteins {
    
    publishDir "${params.newSTRING_rootFolder}", mode: "copy"
    
    
    label "simple_py_process"
    
    conda "/mnt/mnemo5/tao/anaconda3/envs/ipykernel_py3" // seesm configuration here not working , but rather need to be done in configuration file, or need to do both ?
    // debug true //echo true echo directive is depreca
    
    input: 
        val currentSpe_TaxID_ch
        path origProSeqPath
        
    output:
        path "${currentSpe_TaxID_ch}/", type: "dir", emit: currentSpeProSeqPath
        path "${currentSpe_TaxID_ch}ByProteins/", type: "dir", emit: currentSpeProSeqPath_ByProteins
        path "${currentSpe_TaxID_ch}/${currentSpe_TaxID_ch}.fa", type: "file", emit: currentSpe_fastaData
    script:
        
    """
        # cp protein seq of current species  to a new folder and separated them by proteins for later use 
        currentSpeProSeqPath="${currentSpe_TaxID_ch}/" 
        mkdir -p \${currentSpeProSeqPath}
        cp "${origProSeqPath}/${currentSpe_TaxID_ch}.fa" \${currentSpeProSeqPath}  # origProSeqPath, here origProSeqPath is a output channel, the "/" at the end is treated as no, in this case ?
        # create .fai inndex file for samtool faxid later , and  parta fasta files by prpteins,need in phmmer section to 
        currentSpe_fastaData="\${currentSpeProSeqPath}${currentSpe_TaxID_ch}.fa" 
        samtools faidx \${currentSpe_fastaData}
        
        currentSpeProSeqPath_ByProteins="${currentSpe_TaxID_ch}ByProteins/"
        mkdir -p \${currentSpeProSeqPath_ByProteins}
        
        # download file "protein.info.v11.5.txt.gz" for the validation reason later 
        currentSpe_protein_info_filename="${currentSpe_TaxID_ch}.protein.info.v11.5.txt.gz" 
        wget  https://stringdb-static.org/download/protein.info.v11.5/${currentSpe_TaxID_ch}.protein.info.v11.5.txt.gz
        
        
        export PYTHONPATH="${projectDir}/../src/utilities/" # "\${projectDir}/../" not working 
        python ${projectDir}/python_scripts/ParseCurSpeFastaByProteins.py --currentSpe_fastaData \${currentSpe_fastaData} --currentSpeProSeqPath_ByProteins \${currentSpeProSeqPath_ByProteins} --currentSpe_protein_info_filename \${currentSpe_protein_info_filename}

    """
    
}


//remove rudadant proteins for later use 
//removed redundant ones: only the longer sequence was kept  if two sequences were over 95% identical and the alignment covered 90% of the shorter sequence
process prepareSingleMSA_RemoveRedundantProteins {
    
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"  //this folder is just newSTRING_rootFolder
    
    
    label "simple_py_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        //path newSTRING_rootFolder  // when downstreaing process add more content to this folder, this process will be trigered and re-run, how to deal with this problem ?? best maynot seprated new folder and this folder, but this will cause too many folders which I dont like ,!!! ah the sollution is change the publishDir to this folder, but remver this folder from the temporary folder in current process !!!, simpleput , the output chanel of new process dont changle the real output , but only go to the temperory working directory, and the inpuzt channel of new process alway avoid to a big folder . 
        val currentSpe_TaxID_ch
        path currentSpe_fastaData
    output:
        path "${currentSpe_TaxID_ch}withinBlast/", type: "dir", emit: currentSpe_withinBlastPath
        // has to output currentSpeProSeqPath_DB also here, so its acturally moved to publishDir, otherwise its only in current process working directory 
        path "${currentSpe_TaxID_ch}_redundant_proteins.csv", type: "dir", emit: redundant_proteins_csvFile
    script:

    """
        currentSpeProSeqPath_DB="${currentSpe_TaxID_ch}DB/${currentSpe_TaxID_ch}"
        mkdir -p \${currentSpeProSeqPath_DB}
        /mnt/mnemo5/tao/BeeiveProgram/ncbi-blast-2.10.0+/bin/makeblastdb -in ${currentSpe_fastaData} -dbtype "prot" \
        -out \${currentSpeProSeqPath_DB} -parse_seqids
        
        
        currentSpe_withinBlastPath="${currentSpe_TaxID_ch}withinBlast/"
        echo \${currentSpe_withinBlastPath}
        mkdir -p \${currentSpe_withinBlastPath}
        
        /mnt/mnemo5/tao/BeeiveProgram/ncbi-blast-2.10.0+/bin/blastp -num_threads 1 -query ${currentSpe_fastaData} \
         -db \${currentSpeProSeqPath_DB} \
         -out "\${currentSpe_withinBlastPath}all2all.txt" \
         -evalue 1e-6  \
         -outfmt '7 qseqid qaccver  qlen sseqid saccver slen qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovs qcovhsp'
        
        
        redundant_proteins_csvFile="${currentSpe_TaxID_ch}_redundant_proteins.csv"
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/RemoveRedundantProteins.py --redundant_proteins_csvFile \${redundant_proteins_csvFile} --currentSpe_withinBlastPath \${currentSpe_withinBlastPath}
    """
    
}

        
//then preprocess eggNOG othologous group , to make sure for each orthologous group ,only one protein fro one speices 

//this substep take long time , so put it in one single process for easy debugging  
process prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs {
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"
    
    label "large_memory_process"
    
    // debug true //echo true echo directive is deprecated , here too much output, so delete this line 
    
    input: 
        val currentSpe_TaxID_ch
        val current_EggNOG_maxLevel_ch
        path currentSpe_fastaData
        path eggNOG_folder
        path species_file
        path tree_file
    output:
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_orthologs/", emit: currentSpe_currentMaxLevel_orthologs
    
    script: 
    """
        currentSpe_currentMaxLevel_orthologs="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_orthologs/"
        mkdir -p \${currentSpe_currentMaxLevel_orthologs}
        eggNOG_group_folder="${eggNOG_folder}/groups"
        
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/choose_orthologs_STRING11.05.py -o \${currentSpe_currentMaxLevel_orthologs} \
        -i  ${currentSpe_fastaData} -m ${current_EggNOG_maxLevel_ch} \
        -g \${eggNOG_group_folder} \
        -s ${species_file} -t ${tree_file}
    """
}
 

// for each OG group, find their sequence and save them in to one fasta file ,
// first seuquence has to be query protein for downstream filtering 
process prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas {
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"
    
    label "many_cpu_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        val currentSpe_TaxID_ch
        val current_EggNOG_maxLevel_ch
        path currentSpe_currentMaxLevel_orthologs
        path redundant_proteins_csvFile
        path origSTRINGBacteriaProSeqPath
        path currentSpe_fastaData
    output:
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_MiddleData/",type: "dir",  emit: currentSpeMiddleDataPath
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_MiddleData/newsingleMSA_RBH_OrthologousGroup.csv",type: "file", emit: newsingleMSA_RBH_OrthologousGroup_fileName
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_EggNOG_OrthologousGroup_Fa/", type: "dir", emit: currentSpe_OrthologousGroup_Fa_path
    
    script: 
    """
        currentSpeMiddleDataPath="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_MiddleData/"
        mkdir -p \${currentSpeMiddleDataPath} # create the folder to prevent non-existing folder/file problem later
        newsingleMSA_RBH_OrthologousGroup_fileName="\${currentSpeMiddleDataPath}newsingleMSA_RBH_OrthologousGroup.csv"
        
        
 currentSpe_OrthologousGroup_Fa_path="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_EggNOG_OrthologousGroup_Fa/"

       mkdir -p \${currentSpe_OrthologousGroup_Fa_path}
       
       currentSpe_OrthologousGroup_Fa_logpath="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_EggNOG_OrthologousGroup_Fa_log/"
       
       mkdir -p \${currentSpe_OrthologousGroup_Fa_logpath} 
       
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        #first part of this process is fast, so good for debug
        python ${projectDir}/python_scripts/PreprocessEggnogOrthologGroup_collectingOGFastas.py -sfa ${currentSpe_fastaData} \
        -id ${current_EggNOG_maxLevel_ch} -c ${currentSpe_currentMaxLevel_orthologs} \
        -r ${redundant_proteins_csvFile} -f \${newsingleMSA_RBH_OrthologousGroup_fileName} \
        -fa \${currentSpe_OrthologousGroup_Fa_path} -log \${currentSpe_OrthologousGroup_Fa_logpath} \
        -b "${origSTRINGBacteriaProSeqPath}/" -n ${params.small_mp_task_nums} -ut ${params.code_utilities_folder}
        # here add "/" to the currentSpeProSeqPath_ByProteins, as nextflow remove "/" at the end by default
        
    """
}

    


// “SeedAlignment” section , Phmmer choose most similar sequce and hmmalig to do multiple alignment ?
// ???!!! this step seem generate slightly different resutls, because I didnt set random set in phmmer progrom ?
// but seem input currentSpe_OrthologousGroup_Fa_path also slightly differnt ?
// then find that  process chooseOrthologs already give alight different results , this is because when protein have multile og in one spe, current spe randomly choose 1  (level_og_prot have "set" values)?, this is acceptable behaviour as it happens quite rare, ignore this problem for now
process prepareSingleMSA_SeedAlignment {
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"
    
    label "many_cpu_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        val currentSpe_TaxID_ch
        val current_EggNOG_maxLevel_ch
        path currentSpeProSeqPath_ByProteins
        path currentSpe_OrthologousGroup_Fa_path
    output:
         path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_phmmer_results_pident/", emit: currentSpe_phmmer_outPath
        
    script: 
    """        
 currentSpe_phmmer_outPath="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_phmmer_results_pident/"

       mkdir -p \${currentSpe_phmmer_outPath}
       
       currentSpe_phmmer_logPath="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_phmmer_results_pident_log/"
       
       mkdir -p \${currentSpe_phmmer_logPath} 
       
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        # here add "/" to the currentSpeProSeqPath_ByProteins, as nextflow remove "/" at the end by default
        python ${projectDir}/python_scripts/SeedAlignment.py -b "${currentSpeProSeqPath_ByProteins}/" -ofa "${currentSpe_OrthologousGroup_Fa_path}/" \
        -out \${currentSpe_phmmer_outPath} -log \${currentSpe_phmmer_logPath} \
        -n ${params.middle_mp_task_nums} -ut ${params.code_utilities_folder} -ph ${params.phmmer_path}
        
    """
}

     



// # then filter phmmer resutls using
process prepareSingleMSA_SeedAlignment_filtering {
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"
    
    label "large_memory_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        val currentSpe_TaxID_ch
        val current_EggNOG_maxLevel_ch
        path currentSpeProSeqPath_ByProteins
        path origSTRINGBacteriaProSeqPath
        path currentSpe_phmmer_outPath
    output:
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_phmmer_OrthologousGroup_Fa_pident/", emit: currentSpe_phmmer_OrthologousGroup_path
    script: 
    """      
                                            currentSpe_phmmer_OrthologousGroup_path="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_phmmer_OrthologousGroup_Fa_pident/"
        mkdir -p \${currentSpe_phmmer_OrthologousGroup_path}
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/SeedAlignment_filtering.py -b "${currentSpeProSeqPath_ByProteins}/" -pf  \${currentSpe_phmmer_OrthologousGroup_path} \
        -o "${origSTRINGBacteriaProSeqPath}/" -p "${currentSpe_phmmer_outPath}/" -n ${params.small_mp_task_nums}
        
    """
}


//  “ProfileHMM” section, 
// before to make HMM profile from seed alignmt. we first need to perfomacne multiple sequence alignment 
process prepareSingleMSA_ProfileHMM_ClustoMA {
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"
    
    label "manyCPU_largeMemory_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        val currentSpe_TaxID_ch
        val current_EggNOG_maxLevel_ch
        path currentSpe_phmmer_OrthologousGroup_path
    output:
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_phmmer_OrthologousGroup_Fa_pident_ClustoMSA/", emit: currentSpe_ClustoMSA_path
    script: 
    """      

    currentSpe_ClustoMSA_path="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_phmmer_OrthologousGroup_Fa_pident_ClustoMSA/"
        mkdir -p \${currentSpe_ClustoMSA_path}
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/ProfileHMM_ClustoMA.py -pf "${currentSpe_phmmer_OrthologousGroup_path}/" \
        -c \${currentSpe_ClustoMSA_path} -n ${params.large_mp_task_nums} -cp ${params.clustalo_path}
        
    """
}



// now build hmm profile using seed alignmen
process prepareSingleMSA_ProfileHMM_hmmbuild {
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"
    
    label "manyCPU_largeMemory_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        val currentSpe_TaxID_ch
        val current_EggNOG_maxLevel_ch
        path currentSpe_ClustoMSA_path
    output:
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_hmm_profiles_pident/", emit: currentSpe_hmm_profiles_path
        
    script: 
    """      
        currentSpe_hmm_profiles_path="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_hmm_profiles_pident/"
        mkdir -p \${currentSpe_hmm_profiles_path}
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/ProfileHMM_hmmbuild.py -hmm  \${currentSpe_hmm_profiles_path} \
        -c "${currentSpe_ClustoMSA_path}/" -n ${params.large_mp_task_nums} -hmmp ${params.hmmbuild_path}
        
    """
}



// “HMMAllAlignment”  secftion
// now build hmm profile using seed alignmen
process prepareSingleMSA_HMMAllAlignment {
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"
    
    label "manyCPU_largeMemory_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        val currentSpe_TaxID_ch
        val current_EggNOG_maxLevel_ch
        path currentSpe_hmm_profiles_path
        path currentSpe_OrthologousGroup_Fa_path
    output:
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_hmmalign_pident_results/", emit: currentSpe_hmm_align_path
        
    script: 
    """      
        currentSpe_hmm_align_path="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_hmmalign_pident_results/"
        mkdir -p \${currentSpe_hmm_align_path}
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/HMMAllAlignment.py  -hmm "${currentSpe_hmm_profiles_path}/" \
        -a \${currentSpe_hmm_align_path} -ofa "${currentSpe_OrthologousGroup_Fa_path}/" \
        -n ${params.large_mp_task_nums} -align ${params.hmmalign_path}
        
    """
}



// # # Here we need to further remove gaps from single MSA and keep track their position in original seq 
// # here we first remove gaps in ref seq / first protein 
// # then remove seq with many gaps 
// # then remove columns with many gaps 
process prepareSingleMSA_singleMSAGapsFiltering {
    publishDir "${params.newSTRING_rootFolder}",mode: "copy"
    
    label "many_cpu_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        val currentSpe_TaxID_ch
        val current_EggNOG_maxLevel_ch
        path currentSpe_hmm_align_path
        path currentSpe_ClustoMSA_path
        // path currentSpeMiddleDataPath
    output:
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_hmmalign_removeGaps_keepGapPos/", emit: currentSpe_msa_removeGaps_path
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_hmmalign_removeGaps_trackGapPos/", emit: currentSpe_msa_trackGapsPos_path
        path "${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_MSAGapsFiltering/", emit:currentSpeMSAGapsFilteringMetaFolder
        
    script: 
    """      
        currentSpe_msa_removeGaps_path="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_hmmalign_removeGaps_keepGapPos/"
        mkdir -p \${currentSpe_msa_removeGaps_path} 
                currentSpe_msa_trackGapsPos_path="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_newSingleMSA_hmmalign_removeGaps_trackGapPos/"
        mkdir -p \${currentSpe_msa_trackGapsPos_path}
        
        currentSpeMSAGapsFilteringMetaFolder="${currentSpe_TaxID_ch}_EggNOGmaxLevel${current_EggNOG_maxLevel_ch}_MSAGapsFiltering/"
        mkdir -p \${currentSpeMSAGapsFilteringMetaFolder} # this is created to replace old "currentSpeMiddleDataPath" to prevent the current process to change input channel, bad bad, the previous will be trigled when re-run https://github.com/TaoDFang/MNF/issues/103#issuecomment-1694210857
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/singleMSAGapsFiltering.py -a "${currentSpe_hmm_align_path}/" \
        -c "${currentSpe_ClustoMSA_path}/" -tp \${currentSpe_msa_trackGapsPos_path} -rp \${currentSpe_msa_removeGaps_path} \
        -m \${currentSpeMSAGapsFilteringMetaFolder} -n ${params.middle_mp_task_nums}
        
    """
}
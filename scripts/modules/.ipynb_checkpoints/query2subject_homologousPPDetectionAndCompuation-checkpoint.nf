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
    
    // publishDir "${homologous_SeqMappingPath}", mode: "copy" , here the output is aboluste path so we dont need to publish data
    // the reason we do this is that this process need to be run multiple times parallely and need to write to same same folder, 
    // this cause problems 
    // but use abosulut path reduce "explicitely" dependency between process ?, check what is nextflow diagraw looks like here 
    //while we use fromPath to "implicatily" to capture this dependence ?
    // or not use SubjectProSeqPath_ByProtein_ch but a list of pure subject ids ???
    
    
    label "many_cpu_process"
    
    debug true //echo true echo directive is deprecated
    
    
    input: 
        val Query_tuple_ch
        val Subject_tupleList_ch
        // path SubjectSpe_MiddleData_ch
        each QueryProSeqPath_ByProtein_ch
        each SubjectProSeqPath_ByProtein_ch //use each here bacause its the onyl inpuzt channle have multiple element; https://carpentries-incubator.github.io/workflows-nextflow/05-processes-part1/index.html
    // #before using each keyword, the path  offSubjectProSeqPath_ByProteinSubjectProSeqPath_ByProtein: 411476ByProteins/
    // #when use each key workd for , we got full path SubjectProSeqPath_ByProtein: /mnt/mnemo6/tao/nextflow/PPI_Coevolution/STRING_data_11.5/411476ByProteins/, become a value channle now ??
        path homologous_allQuery2SubjectPPIMapping_path

        
//     output:
//         path "${params.homologous_SeqMappingPath}/EggNogMaxLevel2_QuerySpe_ID${params.query_currentSpe_TaxID}andSubjectSpe_ID\${Subject_speID}/", emit: current_homologous_SeqMappingPath
    // Subject_speID is defined with bash script, cause problelm, maybe exract diffrect from by SubjectProSeqPath_ByProtein_ch by grooy command
        
    script:
        
    """
        echo "process domologousPPDetection_SeqMapping: ${SubjectProSeqPath_ByProtein_ch}"
        
        f=\$(basename -- "${SubjectProSeqPath_ByProtein_ch}")  # https://stackoverflow.com/questions/66568781/how-to-call-a-variable-created-in-the-script-in-nextflow
        echo \${f}
        
        [[ "\${f}" =~ ([0-9]+)* ]] #https://www.bashsupport.com/bash/variables/bash/bash_rematch/
            
        Subject_speID=\$(echo \${BASH_REMATCH[1]})
        
        echo "Subject_speID: \${Subject_speID}"
        
    current_homologous_SeqMappingPath="${params.homologous_SeqMappingPath}/EggNogMaxLevel2_QuerySpe_ID${params.query_currentSpe_TaxID}andSubjectSpe_ID\${Subject_speID}/"

        mkdir -p \${current_homologous_SeqMappingPath}
        
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/homologousPPDetection_SeqMapping.py -q ${Query_tuple_ch.join("_")}  \
        -s ${Subject_tupleList_ch.join("_")} \
        -qb "${QueryProSeqPath_ByProtein_ch}/" -sb "${SubjectProSeqPath_ByProtein_ch}/"  -seqM \${current_homologous_SeqMappingPath} \
        -m "${homologous_allQuery2SubjectPPIMapping_path}/" -bp ${params.blastp_path} -n ${params.middle_mp_task_nums}
               
        
    """
    
}


// do protein mapping for homologous pp and single protein mapping, not the one best homologous pp
process homologousPPDetection_allQuery2SubjectPPIMapping_BestHomologous {
    
    publishDir "${params.newSTRING_rootFolder}", mode: "copy"
    
    
    label "simple_py_process"
    
    debug true //echo true echo directive is deprecated
    
    
    input: 
        val Query_tuple_ch
        val Subject_tupleList_ch
        path homologous_SeqMappingPath_ch
        path homologous_allQuery2SubjectPPIMapping_path
        
    output:
        path "${Query_tuple_ch[1]}_EggNOGmaxLevel${Query_tuple_ch[0]}_allQuery2SubjectPPIMapping_singleProteinBlastp/", emit: homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path
        path "${Query_tuple_ch[1]}_EggNOGmaxLevel${Query_tuple_ch[0]}_allQuery2SubjectPPIMapping_bestHomologousPP/" , emit: homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path
    script:
        
    """
        homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path="${Query_tuple_ch[1]}_EggNOGmaxLevel${Query_tuple_ch[0]}_allQuery2SubjectPPIMapping_singleProteinBlastp/" 
        mkdir -p \${homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path}
        
        homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path="${Query_tuple_ch[1]}_EggNOGmaxLevel${Query_tuple_ch[0]}_allQuery2SubjectPPIMapping_bestHomologousPP/" 
        mkdir   -p \${homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path}


        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python   ${projectDir}/python_scripts/homologousPPDetection_allQuery2SubjectPPIMapping_BestHomologous.py  \
        -q ${Query_tuple_ch.join("_")}   -s ${Subject_tupleList_ch.join("_")} \
        -seqM  "${homologous_SeqMappingPath_ch}/"   -m "${homologous_allQuery2SubjectPPIMapping_path}/" \
        -ms \${homologous_allQuery2SubjectPPIMapping_singleProteinBlastp_path} \
        -mb \${homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path}
        
    """
    
}


// // dont includ this step in the computation pipeline but as the post-preporcess steps
// // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_homologousPPDetection_STRINGPhyBalancePhyla.ipynb
// process homologousPPDetection_msa2orig {
        
//     """
//         export PYTHONPATH="${projectDir}/../src/utilities/" 
//         python   ${projectDir}/python_scripts/homologousPPDetection_msa2orig.py  
        
//     """
    
// }


// prepare paired MSA for homologous pp in other species 
process homologousPPDetection_preparePairedMSA {
    
    publishDir "${params.CoEvo_data_folder}", mode: "copy", saveAs: { filename ->  filename.substring(23) } // here the purpose is to only publish subfolder in output path "temp_CoEvo_data_folder/", https://nextflow-io.github.io/patterns/publish-rename-outputs/, https://www.tutorialspoint.com/groovy/groovy_substring.htm
    
    
    label "many_cpu_process"
    
    debug true //echo true echo directive is deprecated
    
    
    input: 
        val Query_tuple_ch
        val Subject_tupleList_ch
        path newSTRING_rootFolder_ch
        path homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path
    output:
        path "temp_CoEvo_data_folder/*", emit: temp_CoEvo_data_folder
       
        
    script:
        
    """
        #create this folder to fit to nextflow logic (the input to next process is better the output of previous process), but the whole folder not just the conent will be move to ? \${params.CoEvo_data_folder}, try with and without /?
        temp_CoEvo_data_folder="temp_CoEvo_data_folder/"
        mkdir -p \${temp_CoEvo_data_folder}

        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/homologousPPDetection_preparePairedMSA.py  \
         -q ${Query_tuple_ch.join("_")}    -s ${Subject_tupleList_ch.join("_")} \
         -f "${newSTRING_rootFolder_ch}/"  -c \${temp_CoEvo_data_folder} \
         -mb "${homologous_allQuery2SubjectPPIMapping_bestHomologousPP_path}/" \
         -nf90  ${params.Nf90_thres}    -n  ${params.large_mp_task_nums}
        
    """
    
}






// !!!!
//for need to do DCA compuation for millions of PPs , if interupt, better not re-start from scratch 
// in this case, better not use output to the temperory working directory and then copy to publishDir
// but direct use final output folder(absolute path) as the input channel or as input parameters 
process homologousPPDetection_ComputeHomologousDCA {
    // publishDir "${params.input_root_folder}",mode: "copy"
    debug true //echo true echo directive is deprecated
    
    
    label "coevolutionComputation_mfDCA_process"
    
    input: 
        val Query_tuple_ch
        val Subject_tupleList_ch
        path newSTRING_rootFolder_ch
        path CoEvo_data_folder_ch
    // output:
    script: 
    """      
        
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 

        python ${projectDir}/python_scripts/homologousPPDetection_ComputeHomologousDCA.py \
        -q ${Query_tuple_ch.join("_")}     -s ${Subject_tupleList_ch.join("_")} \
        -f "${newSTRING_rootFolder_ch}/"   -c "${CoEvo_data_folder_ch}/" \
        -n  ${params.middle_mp_task_nums}
        
    """
}
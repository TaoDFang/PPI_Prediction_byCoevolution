#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// or this one beter ? http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1





include {downLoadOtherRawFiles; downLoadRawFastaFile; prepareFastaDataBySpecies;moveOnlyBacteriaSepcies } from './modules/RawFastaFilesAndMetaData'



workflow {    
    

    println "scriptFile: " + workflow.scriptFile
    println "projectDir: " + workflow.projectDir
    println "launchDir: " + workflow.launchDir
    println "workDir: " + workflow.workDir
    println "configFiles: " + workflow.configFiles
    
    
    // ************Download all metadata and sequencing data**********
    
    downLoadOtherRawFiles_ch=downLoadOtherRawFiles()
    
    //downLoadRawFastaFile_ch=downLoadRawFastaFile(RawData_Folder_ch)
    downLoadRawFastaFile_ch=downLoadRawFastaFile()
    //downLoadRawFastaFile.out.view()
    // downLoadRawFastaFile_ch.rawFasta_file.view()
  
    prepareFastaDataBySpecies_ch=prepareFastaDataBySpecies(downLoadRawFastaFile_ch.rawFasta_file)
    prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder.view()
    //prepareFastaDataBySpecies.out.view()
    
    // here use channel not params.STRING_fastaBySpecies_Folder directly to triger next process ""
    //STRING_fastaBySpecies_Folder_ch=Channel.fromPath(params.STRING_fastaBySpecies_Folder,type:'dir')
    moveOnlyBacteriaSepcies_ch=moveOnlyBacteriaSepcies(downLoadOtherRawFiles_ch.species_file,prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder)
    //println "moveOnlyBacteriaSepcies.out.view: " + moveOnlyBacteriaSepcies.out.view()
    
    
    
    // ************Prepare single MSA **********
    //http://localhost:8206/lab/workspaces/auto-I/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
    //http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.py
    // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
    
    // ************Prepare paired MSA **********
    prepareSingleMSA_ParseCurSpeFastaByProteins_ch=prepareSingleMSA_ParseCurSpeFastaByProteins(prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder)
    prepareSingleMSA_ParseCurSpeFastaByProteins_ch.newSTRING_rootFolder.view()
    
    prepareSingleMSA_RemoveRedundantProteins_ch=prepareSingleMSA_RemoveRedundantProteins(prepareSingleMSA_ParseCurSpeFastaByProteins_ch.newSTRING_rootFolder,prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpe_fastaData)
    
    
    prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs_ch=prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs(prepareSingleMSA_ParseCurSpeFastaByProteins_ch.newSTRING_rootFolder,prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpe_fastaData,downLoadOtherRawFiles_ch.eggNOG_folder,downLoadOtherRawFiles_ch.species_file,downLoadOtherRawFiles_ch.species_tree_file)
    
    prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch=prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas(prepareSingleMSA_ParseCurSpeFastaByProteins_ch.newSTRING_rootFolder,prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs_ch.currentSpe_currentMaxLevel_orthologs,prepareSingleMSA_RemoveRedundantProteins_ch.redundant_proteins_csvFile,moveOnlyBacteriaSepcies_ch.STRING_fastaByBacteriaSpecies_Folder)
    // ************Compute DCA  **********
    
    
}



// ************Prepare single MSA **********

//section ParseEcoliFastaByProtein, need in phmmer section to 
process prepareSingleMSA_ParseCurSpeFastaByProteins {
    
    publishDir "${params.PPI_Coevolution}", mode: "copy"
    
    
    label "simple_py_process"
    
    conda "/mnt/mnemo5/tao/anaconda3/envs/ipykernel_py3" // seesm configuration here not working , but rather need to be done in configuration file, or need to do both ?
    // debug true //echo true echo directive is depreca
    
    input: 
        path origProSeqPath
        
    output:
        path "STRING_data_11.5/",type: "dir", emit: newSTRING_rootFolder
        path "STRING_data_11.5/${params.currentSpe_TaxID}/", type: "dir", emit: currentSpeProSeqPath
        path "STRING_data_11.5/${params.currentSpe_TaxID}ByProteins/", type: "dir", emit: currentSpeProSeqPath_ByProteins
        path "STRING_data_11.5/${params.currentSpe_TaxID}/${params.currentSpe_TaxID}.fa", type: "file", emit: currentSpe_fastaData
    script:
        
        //newSTRING_rootFolder="${params.PPI_Coevolution}/STRING_data_11.5/" # here do not use params.PPI_Coevolution to avoid abolute path 
    """
        # cp protein seq of current species  to a new folder and separated them by proteins for later use 
        newSTRING_rootFolder="STRING_data_11.5/" #define same in the output to export folder for downstreaming process
        mkdir -p \${newSTRING_rootFolder}
        currentSpeProSeqPath="\${newSTRING_rootFolder}${params.currentSpe_TaxID}/" 
        mkdir -p \${currentSpeProSeqPath}
        cp "${origProSeqPath}/${params.currentSpe_TaxID}.fa" \${currentSpeProSeqPath}  # origProSeqPath, here origProSeqPath is a output channel, the "/" at the end is treated as no, in this case ?
        # create .fai inndex file for samtool faxid later , and  parta fasta files by prpteins,need in phmmer section to 
        currentSpe_fastaData="\${currentSpeProSeqPath}${params.currentSpe_TaxID}.fa" 
        samtools faidx \${currentSpe_fastaData}
        
        currentSpeProSeqPath_ByProteins="\${newSTRING_rootFolder}${params.currentSpe_TaxID}ByProteins/"
        mkdir -p \${currentSpeProSeqPath_ByProteins}
        
        # download file "protein.info.v11.5.txt.gz" for the validation reason later 
        currentSpe_protein_info_filename="\${newSTRING_rootFolder}${params.currentSpe_TaxID}.protein.info.v11.5.txt.gz" 
        wget  https://stringdb-static.org/download/protein.info.v11.5/${params.currentSpe_TaxID}.protein.info.v11.5.txt.gz -P \${newSTRING_rootFolder}
        
        
        export PYTHONPATH="${projectDir}/../src/utilities/" # "\${projectDir}/../" not working 
        python ${projectDir}/python_scripts/ParseCurSpeFastaByProteins.py --currentSpe_fastaData \${currentSpe_fastaData} --currentSpeProSeqPath_ByProteins \${currentSpeProSeqPath_ByProteins} --currentSpe_protein_info_filename \${currentSpe_protein_info_filename}

    """
    
}


//remove rudadant proteins for later use 
//removed redundant ones: only the longer sequence was kept  if two sequences were over 95% identical and the alignment covered 90% of the shorter sequence
process prepareSingleMSA_RemoveRedundantProteins {
    
    publishDir "${params.PPI_Coevolution}",mode: "copy"
    
    
    label "simple_py_process"
    
    // debug true //echo true echo directive is deprecated
    
    input: 
        path newSTRING_rootFolder  // when downstreaing process add more content to this folder, this process will be trigered and re-run, how to deal with this problem ?? best maynot seprated new folder and this folder, but this will cause too many folders which I dont like ,!!! ah the sollution is change the publishDir to this folder, but remver this folder from the temporary folder in current process !!!
        path currentSpe_fastaData
    output:
        path "${newSTRING_rootFolder}/${params.currentSpe_TaxID}withinBlast/", type: "dir", emit: currentSpe_withinBlastPath
        // has to output currentSpeProSeqPath_DB also here, so its acturally moved to publishDir, otherwise its only in current process working directory 
        path "${newSTRING_rootFolder}/${params.currentSpe_TaxID}_redundant_proteins.csv", type: "dir", emit: redundant_proteins_csvFile
    script:

    """
        currentSpeProSeqPath_DB="${newSTRING_rootFolder}/${params.currentSpe_TaxID}DB/${params.currentSpe_TaxID}"
        mkdir -p \${currentSpeProSeqPath_DB}
        /mnt/mnemo5/tao/BeeiveProgram/ncbi-blast-2.10.0+/bin/makeblastdb -in ${currentSpe_fastaData} -dbtype "prot" \
        -out \${currentSpeProSeqPath_DB} -parse_seqids
        
        
        currentSpe_withinBlastPath="${newSTRING_rootFolder}/${params.currentSpe_TaxID}withinBlast/"
        echo \${currentSpe_withinBlastPath}
        mkdir -p \${currentSpe_withinBlastPath}
        
        /mnt/mnemo5/tao/BeeiveProgram/ncbi-blast-2.10.0+/bin/blastp -num_threads 1 -query ${currentSpe_fastaData} \
         -db \${currentSpeProSeqPath_DB} \
         -out "\${currentSpe_withinBlastPath}all2all.txt" \
         -evalue 1e-6  \
         -outfmt '7 qseqid qaccver  qlen sseqid saccver slen qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovs qcovhsp'
        
        
        redundant_proteins_csvFile="${newSTRING_rootFolder}/${params.currentSpe_TaxID}_redundant_proteins.csv"
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/RemoveRedundantProteins.py --redundant_proteins_csvFile \${redundant_proteins_csvFile} --currentSpe_withinBlastPath \${currentSpe_withinBlastPath}
    """
    
}

        
//then preprocess eggNOG othologous group , to make sure for each orthologous group ,only one protein fro one speices 

//this substep take long time , so put it in one single process for easy debugging  
process prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs {
    publishDir "${params.PPI_Coevolution}",mode: "copy"
    
    label "large_memory_process"
    
    // debug true //echo true echo directive is deprecated , here too much output, so delete this line 
    
    input: 
        path newSTRING_rootFolder
        path currentSpe_fastaData
        path eggNOG_folder
        path species_file
        path tree_file
    output:
        path "${newSTRING_rootFolder}/${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_orthologs/", emit: currentSpe_currentMaxLevel_orthologs
    
    script: 
    """
        currentSpe_currentMaxLevel_orthologs="${newSTRING_rootFolder}/${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_orthologs/"
        mkdir -p \${currentSpe_currentMaxLevel_orthologs}
        eggNOG_group_folder="${eggNOG_folder}/groups"
        
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        python ${projectDir}/python_scripts/choose_orthologs_STRING11.05.py -o \${currentSpe_currentMaxLevel_orthologs} \
        -i  ${currentSpe_fastaData} -m ${params.current_EggNOG_maxLevel} \
        -g \${eggNOG_group_folder} \
        -s ${species_file} -t ${tree_file}
    """
}
 

// for each OG group, find their sequence and save them in to one fasta file ,
// first seuquence has to be query protein for downstream filtering 
process prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas {
    publishDir "${params.PPI_Coevolution}",mode: "copy"
    
    label "many_cpu_process"
    
    debug true //echo true echo directive is deprecated
    
    input: 
        path newSTRING_rootFolder
        path currentSpe_currentMaxLevel_orthologs
        path redundant_proteins_csvFile
        path origSTRINGBacteriaProSeqPath
    output:
        path "${newSTRING_rootFolder}/${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_MiddleData/",type: "dir",  emit: currentSpeMiddleDataPath
        path "${newSTRING_rootFolder}/${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_MiddleData/newsingleMSA_RBH_OrthologousGroup.csv",type: "file", emit: newsingleMSA_RBH_OrthologousGroup_fileName
        path "${newSTRING_rootFolder}/${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_newSingleMSA_EggNOG_OrthologousGroup_Fa/", type: "dir", emit: currentSpe_OrthologousGroup_Fa_path
    
    script: 
    """
        currentSpeMiddleDataPath="${newSTRING_rootFolder}/${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_MiddleData/"
        mkdir -p \${currentSpeMiddleDataPath} # create the folder to prevent non-existing folder/file problem later
        newsingleMSA_RBH_OrthologousGroup_fileName="\${currentSpeMiddleDataPath}newsingleMSA_RBH_OrthologousGroup.csv"
        
        
 currentSpe_OrthologousGroup_Fa_path="${newSTRING_rootFolder}/${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_newSingleMSA_EggNOG_OrthologousGroup_Fa/"

       mkdir -p \${currentSpe_OrthologousGroup_Fa_path}
       
       currentSpe_OrthologousGroup_Fa_logpath="${newSTRING_rootFolder}/${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_newSingleMSA_EggNOG_OrthologousGroup_Fa_log/"
       
       mkdir -p \${currentSpe_OrthologousGroup_Fa_logpath} 
       
        
        export PYTHONPATH="${projectDir}/../src/utilities/" 
        #first part of this process is fast, so good for debug
        python ${projectDir}/python_scripts/PreprocessEggnogOrthologGroup_collectingOGFastas.py -sf ${newSTRING_rootFolder} \
        -id ${params.current_EggNOG_maxLevel} -c ${currentSpe_currentMaxLevel_orthologs} \
        -r ${redundant_proteins_csvFile} -f \${newsingleMSA_RBH_OrthologousGroup_fileName} \
        -fa \${currentSpe_OrthologousGroup_Fa_path} -log \${currentSpe_OrthologousGroup_Fa_logpath} \
        -b ${origSTRINGBacteriaProSeqPath} -n ${params.small_mp_task_nums} -ut ${params.code_utilities_folder}
        
        
    """
}

    


/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run PairedMSA_preprocessing_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

*/

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// or this one beter ? http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1





include { downLoadRawFastaFile; prepareFastaDataBySpecies;moveOnlyBacteriaSepcies } from './modules/RawFastaFilesAndMetaData'



workflow {    
    // ************Download all metadata and sequencing data**********
    
    //downLoadRawFastaFile_ch=downLoadRawFastaFile(RawData_Folder_ch)
    downLoadRawFastaFile_ch=downLoadRawFastaFile()
    //downLoadRawFastaFile.out.view()
    // downLoadRawFastaFile_ch.rawFasta_file.view()
  
    
    prepareFastaDataBySpecies_ch=prepareFastaDataBySpecies(downLoadRawFastaFile_ch.rawFasta_file)
    //prepareFastaDataBySpecies.out.view()
    
    // here use channel not params.STRING_fastaBySpecies_Folder directly to triger next process ""
    //STRING_fastaBySpecies_Folder_ch=Channel.fromPath(params.STRING_fastaBySpecies_Folder,type:'dir')
    moveOnlyBacteriaSepcies(downLoadRawFastaFile_ch.STR_species_mem,prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder)
    moveOnlyBacteriaSepcies.out.view()
    
    
    // ************Prepare single MSA **********
    //http://localhost:8206/lab/workspaces/auto-I/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
    //http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.py
    // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
    
    // ************Prepare paired MSA **********
    prepareSingleMSA_ParseCurSpeFastaByProteins_ch=prepareSingleMSA_ParseCurSpeFastaByProteins(prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder)
    prepareSingleMSA_ParseCurSpeFastaByProteins_ch.view()
    
    
    
    
    // ************Compute DCA  **********
    
    
    
}





// ************Prepare single MSA **********

//here prepare folder and raw files necesseary for the python code 
process prepareSingleMSA_ParseCurSpeFastaByProteins {
    publishDir "${params.PPI_Coevolution}", mode: "copy"
    
    conda "/mnt/mnemo5/tao/anaconda3/envs/ipykernel_py3" 
    debug true //echo true echo directive is depreca
    
    label "simple_process"
    
    input: 
        path: origProSeqPath
        
    output:
        path: "STRING_data_11.5/",type: "dir", emit: newSTRING_rootFolder
        path: "STRING_data_11.5/${params.currentSpe_TaxID}/", type: "dir", emit: currentSpeProSeqPath
        path: "STRING_data_11.5/${params.currentSpe_TaxID}ByProteins/", type: "dir", emit: currentSpeProSeqPath_ByProteins
    script:
        
        //newSTRING_rootFolder="${params.PPI_Coevolution}/STRING_data_11.5/" # here do not use params.PPI_Coevolution to avoid abolute path 
    """
        # cp protein seq of current species  to a new folder and separated them by proteins for later use 
        newSTRING_rootFolder="STRING_data_11.5/" #define same in the output to export folder for downstreaming process
        mkdir -p \${newSTRING_rootFolder}
        currentSpeProSeqPath="\${newSTRING_rootFolder}/${params.currentSpe_TaxID}/" 
        mkdir -p \${currentSpeProSeqPath}
        cp "${origProSeqPath}/${params.currentSpe_TaxID}.fa}" \${currentSpeProSeqPath}
        # create .fai inndex file for samtool faxid later , and  parta fasta files by prpteins,need in phmmer section to 
        currentSpe_fastaData="\${newSTRING_rootFolder}/${params.currentSpe_TaxID}.fa" 
        samtools faidx \${currentSpe_fastaData}
        
        currentSpeProSeqPath_ByProteins="\${newSTRING_rootFolder}/${params.currentSpe_TaxID}ByProteins/"
        mkdir -p \${currentSpeProSeqPath_ByProteins}
        
        # download file "protein.info.v11.5.txt.gz" for the validation reason later 
        currentSpe_protein_info_filename="\${newSTRING_rootFolder}/${params.currentSpe_TaxID}.protein.info.v11.5.txt.gz" 
        mkdir -p \${currentSpe_protein_info_filename}
        wget  https://stringdb-static.org/download/protein.info.v11.5/${params.currentSpe_TaxID}.protein.info.v11.5.txt.gz -P \${newSTRING_rootFolder}
        
        python ${projectDir}/python_scripts/prepareFastaDataBySpecies.py --currentSpe_fastaData \${currentSpe_fastaData} --currentSpeProSeqPath_ByProteins \${currentSpeProSeqPath_ByProteins} --currentSpe_protein_info_filename \${currentSpe_protein_info_filename}
        



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

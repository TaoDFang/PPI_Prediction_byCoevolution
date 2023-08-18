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
    
    
    
    // ************Prepare paired MSA **********
    
    
    
    // ************Compute DCA  **********
    
    
    
}





// ************Prepare single MSA **********







/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run PairedMSA_preprocessing_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

*/

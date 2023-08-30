#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// or this one beter ? http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1




include {downLoadOtherRawFiles; downLoadRawFastaFile; prepareFastaDataBySpecies;moveOnlyBacteriaSepcies } from './modules/RawFastaFilesAndMetaData.nf' //when this folder is in /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts



workflow RawFastaFilesAndMetaData_workflow{    
    

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
    
    
    
}


    
/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training

when this folder is in /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run RawFastaFilesAndMetaData_workflow.nf -entry RawFastaFilesAndMetaData_workflow -params-file wc-params.json -c nextflow.config -resume


with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

this pipeline only run once irregular of species and phylum 
*/


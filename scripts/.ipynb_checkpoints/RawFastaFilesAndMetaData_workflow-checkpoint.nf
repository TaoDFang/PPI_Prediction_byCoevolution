#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// or this one beter ? http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1




include {test_configuration;downLoadOtherRawFiles; downLoadRawFastaFile; prepareFastaDataBySpecies;moveOnlyBacteriaSepcies } from './modules/RawFastaFilesAndMetaData.nf' //when this folder is in /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts



workflow RawFastaFilesAndMetaData_workflow{    
    main: 
        println "scriptFile: " + workflow.scriptFile
        println "projectDir: " + workflow.projectDir
        println "launchDir: " + workflow.launchDir
        println "workDir: " + workflow.workDir
        println "configFiles: " + workflow.configFiles


        // ************Download all metadata and sequencing data**********

        test_configuration()
    
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
    emit: 
        species_file = downLoadOtherRawFiles_ch.species_file
        species_tree_file = downLoadOtherRawFiles_ch.species_tree_file
        eggNOG_folder = downLoadOtherRawFiles_ch.eggNOG_folder
        rawFasta_file = downLoadRawFastaFile_ch.rawFasta_file
        STRING_fastaBySpecies_Folder = prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder
        STRING_fastaByBacteriaSpecies_Folder = moveOnlyBacteriaSepcies_ch.STRING_fastaByBacteriaSpecies_Folder

}



    
/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training

when this folder is in /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run RawFastaFilesAndMetaData_workflow.nf -entry RawFastaFilesAndMetaData_workflow -params-file wc-params.json -c nextflow.config -resume

one hpc cluster 
cd /home/tfang/PPI_Prediction_byCoevolution/scripts
conda activate /data/tfang/conda-envs/nf-training
then run from login ndoe (not from interactive session)
nextflow run RawFastaFilesAndMetaData_workflow.nf -entry RawFastaFilesAndMetaData_workflow  -c nextflow.config -profile slurm  -resume

new local configuration file  with singularity 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run RawFastaFilesAndMetaData_workflow.nf  -entry RawFastaFilesAndMetaData_workflow -c nextflow.config -profile singularity  -with-singularity /mnt/mnemo5/tao/singularity_containers/PPICoe.sif -resume


with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

this pipeline only run once irregular of species and phylum 
*/


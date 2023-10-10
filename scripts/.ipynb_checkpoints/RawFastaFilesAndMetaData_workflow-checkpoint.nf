#!/usr/bin/env nextflow

nextflow.enable.dsl=2



include {test_configuration;downLoadOtherRawFiles; downLoadRawFastaFile; prepareFastaDataBySpecies;moveOnlyBacteriaSepcies } from './modules/RawFastaFilesAndMetaData.nf' 


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

        downLoadRawFastaFile_ch=downLoadRawFastaFile()

        prepareFastaDataBySpecies_ch=prepareFastaDataBySpecies(downLoadRawFastaFile_ch.rawFasta_file)
        prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder.view()

        moveOnlyBacteriaSepcies_ch=moveOnlyBacteriaSepcies(downLoadOtherRawFiles_ch.species_file,prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder)

    emit: 
        species_file = downLoadOtherRawFiles_ch.species_file
        species_tree_file = downLoadOtherRawFiles_ch.species_tree_file
        eggNOG_folder = downLoadOtherRawFiles_ch.eggNOG_folder
        rawFasta_file = downLoadRawFastaFile_ch.rawFasta_file
        STRING_fastaBySpecies_Folder = prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder
        STRING_fastaByBacteriaSpecies_Folder = moveOnlyBacteriaSepcies_ch.STRING_fastaByBacteriaSpecies_Folder

}



    
/*

conda activate nf-training


cd ~/PPI_Prediction_byCoevolution/scripts
nextflow run RawFastaFilesAndMetaData_workflow.nf -entry RawFastaFilesAndMetaData_workflow -params-file wc-params.json -c nextflow.config -resume



new local configuration file  with singularity 
conda activate nf-training
cd ~/PPI_Prediction_byCoevolution/scripts
nextflow run RawFastaFilesAndMetaData_workflow.nf  -entry RawFastaFilesAndMetaData_workflow -c nextflow.config -profile singularity  -with-singularity /mnt/mnemo5/tao/singularity_containers/PPICoe.sif -resume


with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

this pipeline only run once irregular of species and phylum 
*/


#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.CoEvo_data_folder="${params.PPI_Coevolution}/CoEvo_data_STRING11.5/"
params.input_root_folder="${params.CoEvo_data_folder}allPPI_${params.query_currentSpe_TaxID}_EggNOGmaxLevel${params.query_current_EggNOG_maxLevel}_eggNOGfilteredData"// in the next 


include {RawFastaFilesAndMetaData_workflow } from './RawFastaFilesAndMetaData_workflow.nf'
include {prepareSingleMSA_workflow} from './prepareSingleMSA_workflow.nf'


include {preparePairedMSA_oneRunWithNf90;preparePairedMSA_removeHomologousPairs} from  "./modules/preparePairedMSA.nf"

workflow Query_PairedMSA_preprocessing_workflow{    
    take:
        currentSpe_TaxID_ch
        currentSpeMSAGapsFilteringMetaFolder
        currentSpe_msa_removeGaps_path
    main:

        println "scriptFile: " + workflow.scriptFile
        println "projectDir: " + workflow.projectDir
        println "launchDir: " + workflow.launchDir
        println "workDir: " + workflow.workDir
        println "configFiles: " + workflow.configFiles




        // ************Prepare paired MSA **********
        preparePairedMSA_oneRunWithNf90_ch=preparePairedMSA_oneRunWithNf90(currentSpe_TaxID_ch,currentSpeMSAGapsFilteringMetaFolder,currentSpe_msa_removeGaps_path)
        preparePairedMSA_removeHomologousPairs_ch=preparePairedMSA_removeHomologousPairs(preparePairedMSA_oneRunWithNf90_ch.pairedMSA_Nf90_csv,preparePairedMSA_oneRunWithNf90_ch.pairedMSA_Nf90_folder)
    
    emit:
        
        pairedMSA_Nf90_folder=preparePairedMSA_oneRunWithNf90_ch.pairedMSA_Nf90_folder
        PPIInfoBeforeCoEvoComp_csv=preparePairedMSA_removeHomologousPairs_ch.PPIInfoBeforeCoEvoComp_csv
    
}


workflow {
    
    RawFastaFilesAndMetaData_workflow() 
    
    query_prepareSingleMSA_workflow=prepareSingleMSA_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                         Channel.value("${params.query_current_EggNOG_maxLevel}"),
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                              RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                              RawFastaFilesAndMetaData_workflow.out.species_file,
                              RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                             )
    
    Query_PairedMSA_preprocessing_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                                           query_prepareSingleMSA_workflow.currentSpeMSAGapsFilteringMetaFolder,
                                           query_prepareSingleMSA_workflow.currentSpe_msa_removeGaps_path)
}


/*

conda activate nf-training
cd ~/PPI_Prediction_byCoevolution/scripts
nextflow run Query_PairedMSA_preprocessing_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"


*/


#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.CoEvo_data_folder="${params.PPI_Coevolution}/CoEvo_data_STRING11.5/"
params.input_root_folder="${params.CoEvo_data_folder}allPPI_${params.query_currentSpe_TaxID}_EggNOGmaxLevel${params.query_current_EggNOG_maxLevel}_eggNOGfilteredData"// in the 
params.DCA_coevolutoin_path="${params.input_root_folder}/coevolutoin_result_DCA/"



include {RawFastaFilesAndMetaData_workflow } from './RawFastaFilesAndMetaData_workflow.nf'
include {prepareSingleMSA_workflow} from './prepareSingleMSA_workflow.nf'
include {Query_PairedMSA_preprocessing_workflow} from  "./Query_PairedMSA_preprocessing_workflow.nf"

include {coevolutionComputation_mfDCA_createdFolder;coevolutionComputation_mfDCA_preparaIndexFile;coevolutionComputation_mfDCA_parallel} from "./modules/coevolutionComputation"




workflow Query_coevolutionComputation_workflow{    
    take:
        
        DCA_coevolutoin_path_ch
        IndexDCA_idxCH
        DCA_blockNum_ch
        currentSpeMSAGapsFilteringMetaFolder
        PPIInfoBeforeCoEvoComp_csv
        pairedMSA_Nf90_folder

    main:
        
        println "scriptFile: " + workflow.scriptFile
        println "projectDir: " + workflow.projectDir
        println "launchDir: " + workflow.launchDir
        println "workDir: " + workflow.workDir
        println "configFiles: " + workflow.configFiles


        // ************Compute DCA  **********
        // !!!!
        //we need to do DCA compuation for millions of PPs , if interupt, better not re-start from scratch 
        // in this case, better not use output to the temperory working directory and then copy to publishDir
        // but direct use final output folder as the input channel or as input parameters 
        // DCA_coevolutoin_path_ch=Channel.fromPath(params.DCA_coevolutoin_path,type:'dir')
        coevolutionComputation_mfDCA_createdFolder()
        coevolutionComputation_mfDCA_preparaIndexFile_ch=coevolutionComputation_mfDCA_preparaIndexFile(DCA_coevolutoin_path_ch,currentSpeMSAGapsFilteringMetaFolder,PPIInfoBeforeCoEvoComp_csv,pairedMSA_Nf90_folder,DCA_blockNum_ch)
        coevolutionComputation_mfDCA_parallel_ch=coevolutionComputation_mfDCA_parallel(coevolutionComputation_mfDCA_preparaIndexFile_ch.IndexDCA_coevolutoin_path,IndexDCA_idxCH)

}






workflow {

    // ************Download all metadata and sequencing data**********
    RawFastaFilesAndMetaData_workflow() 
    
   // ************Prepare single MSA **********
    query_prepareSingleMSA_workflow=prepareSingleMSA_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                         Channel.value("${params.query_current_EggNOG_maxLevel}"),
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                              RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                              RawFastaFilesAndMetaData_workflow.out.species_file,
                              RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                             )

    // ************Prepare paired MSA **********
    Query_PairedMSA_preprocessing_workflow_ch=Query_PairedMSA_preprocessing_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                                           query_prepareSingleMSA_workflow.currentSpeMSAGapsFilteringMetaFolder,
                                           query_prepareSingleMSA_workflow.currentSpe_msa_removeGaps_path)
    
    
    // ************Compute DCA  **********
    // !!!!
    //we need to do DCA compuation for millions of PPs , if interupt, better not re-start from scratch 
    // in this case, better not use output to the temperory working directory and then copy to publishDir
    // but direct use final output folder as the input channel or as input parameters 
    // DCA_coevolutoin_path_ch=Channel.fromPath(params.DCA_coevolutoin_path,type:'dir')
    DCA_coevolutoin_path_ch=Channel.fromPath(params.DCA_coevolutoin_path,type:'dir')
    IndexDCA_idxCH=Channel.of( 0..params.DCA_blockNum-1 )  // since in the python index starrt from 0
    IndexDCA_idxCH.view()
    DCA_blockNum_ch=Channel.value( "${params.DCA_blockNum}" )
   Query_coevolutionComputation_workflow(DCA_coevolutoin_path_ch,
                                         IndexDCA_idxCH,
                                         DCA_blockNum_ch,
                                        query_prepareSingleMSA_workflow.currentSpeMSAGapsFilteringMetaFolder,
                                         Query_PairedMSA_preprocessing_workflow_ch.PPIInfoBeforeCoEvoComp_csv,
                                         Query_PairedMSA_preprocessing_workflow_ch.pairedMSA_Nf90_folder,
                                        )
    
    

}



    

/*

conda activate nf-training
cd ~/PPI_Prediction_byCoevolution/scripts
nextflow run Query_coevolutionComputation_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

new local configuration file 
conda activate nf-training
cd ~/PPI_Prediction_byCoevolution/scripts
nextflow run Query_coevolutionComputation_workflow.nf  -c nextflow.config -profile standard  -resume

new local configuration file  with singularity 
conda activate nf-training
cd ~/PPI_Prediction_byCoevolution/scripts
nextflow run Query_coevolutionComputation_workflow.nf  -c nextflow.config -profile singularity   -resume


*/

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// or this one beter ? http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1


// add this param here after inlcude as params.input_root_folder is from Query_PairedMSA_preprocessing_workflow.nf
// not working , need to add directly here params.input_root_folder 
// this will overwrite the definition in Query_PairedMSA_preprocessing_workflow.nf,
// we still have definition of params.input_root_folder in Query_PairedMSA_preprocessing_workflow.nf just incase we need to run it along 
params.CoEvo_data_folder="${params.PPI_Coevolution}/CoEvo_data_STRING11.5/"
params.input_root_folder="${params.CoEvo_data_folder}allPPI_${params.query_currentSpe_TaxID}_EggNOGmaxLevel${params.query_current_EggNOG_maxLevel}_eggNOGfilteredData"// in the next nextflow pipeline, by default its query species folder 
params.DCA_coevolutoin_path="${params.input_root_folder}/coevolutoin_result_DCA/"
// params.IndexDCA_coevolutoin_path="${params.input_root_folder}/coevolutoin_computation_IndexDCA/"


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

            // coevolutionComputation_mfDCA_ch=coevolutionComputation_mfDCA(DCA_coevolutoin_path_ch,currentSpeMSAGapsFilteringMetaFolder,PPIInfoBeforeCoEvoComp_csv,pairedMSA_Nf90_folder)

    
}






workflow {

    // ************Download all metadata and sequencing data**********
    RawFastaFilesAndMetaData_workflow() 
    
   // ************Prepare single MSA **********
    //http://localhost:8206/lab/workspaces/auto-I/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
    //http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.py
    query_prepareSingleMSA_workflow=prepareSingleMSA_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                         Channel.value("${params.query_current_EggNOG_maxLevel}"),
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                              RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                              RawFastaFilesAndMetaData_workflow.out.species_file,
                              RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                              RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                             )

    // ************Prepare paired MSA **********
    // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
    // here for the because the output generate millions of output,  to prevent completely re-run of progrom in case of accidentaly interuption?
    // some input channle should be the final output path, not the the temperaoy woring directory of previous path 
    // but actully if last step finished correctly , its okay, more the problem of last step compute DCA
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
    // IndexDCA_coevolutoin_path_ch=Channel.fromPath(params.IndexDCA_coevolutoin_path,type:'dir')
   Query_coevolutionComputation_workflow(DCA_coevolutoin_path_ch,
                                         IndexDCA_idxCH,
                                         DCA_blockNum_ch,
                                        query_prepareSingleMSA_workflow.currentSpeMSAGapsFilteringMetaFolder,
                                         Query_PairedMSA_preprocessing_workflow_ch.PPIInfoBeforeCoEvoComp_csv,
                                         Query_PairedMSA_preprocessing_workflow_ch.pairedMSA_Nf90_folder,
                                        )
    
    
    // here notice, when the workflow/process is from left of "=", we use .currentSpeMSAGapsFilteringMetaFolder
    // but if without assignment, we use .out.currentSpeMSAGapsFilteringMetaFolder,

}



    
// mkdir -p ${DCA_coevolutoin_path}
// mkdir -p ${params.DCA_coevolutoin_path}


// then homologous DCAs??
// one start fiel to check http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/Compute_allPPI_homologousPPDetection.ipynb

    

//         DCA_coevolutoin_path="coevolutoin_result_DCA/"
//         MI_coevolutoin_path="coevolutoin_result_MI/"

//         mkdir -p \${DCA_coevolutoin_path}
//         mkdir -p \${MI_coevolutoin_path}
/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run Query_coevolutionComputation_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

new local configuration file 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run Query_coevolutionComputation_workflow.nf  -c nextflow.config -profile standard  -resume

new local configuration file  with singularity 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run Query_coevolutionComputation_workflow.nf  -c nextflow.config -profile singularity   -resume



one hpc cluster 
cd /home/tfang/PPI_Prediction_byCoevolution/scripts
conda activate /data/tfang/conda-envs/nf-training
then run from login ndoe (now from srun interactive session)
nextflow run Query_coevolutionComputation_workflow.nf  -c nextflow.config -profile slurm  -resume


*/

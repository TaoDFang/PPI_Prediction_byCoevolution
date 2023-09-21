#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// or this one beter ? http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

// https://github.com/TaoDFang/MNF/issues/98
// https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1

//something changed in vscode, check data change in jupyter

params.CoEvo_data_folder="${params.PPI_Coevolution}/CoEvo_data_STRING11.5/"
params.input_root_folder="${params.CoEvo_data_folder}allPPI_${params.query_currentSpe_TaxID}_EggNOGmaxLevel${params.query_current_EggNOG_maxLevel}_eggNOGfilteredData"// in the next nextflow pipeline, by default its query species folder 

//params.DCA_coevolutoin_path="${params.input_root_folder}/coevolutoin_result_DCA/"


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
        // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
        // here for the because the output generate millions of output,  to prevent completely re-run of progrom in case of accidentaly interuption?
        // some input channle should be the final output path, not the the temperaoy woring directory of previous path 
        // but actully if last step finished correctly , its okay, more the problem of last step compute DCA
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
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run Query_PairedMSA_preprocessing_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

to run for other species and phylum
nextflow run Query_PairedMSA_preprocessing_workflow.nf --currentSpe_TaxID "1274374" --current_EggNOG_maxLevel "1239" -params-file wc-params.json -c nextflow.config -resume


one hpc cluster 
cd /home/tfang/PPI_Prediction_byCoevolution/scripts
conda activate /data/tfang/conda-envs/nf-training
then run from login ndoe (now from srun interactive session)
nextflow run Query_PairedMSA_preprocessing_workflow.nf  -c nextflow.config -profile slurm  -resume


*/


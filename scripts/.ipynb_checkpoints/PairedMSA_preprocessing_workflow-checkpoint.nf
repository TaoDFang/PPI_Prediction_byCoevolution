#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// or this one beter ? http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1


params.newSTRING_rootFolder="${params.PPI_Coevolution}/STRING_data_11.5" //this folder is just newSTRING_rootFolder
params.CoEvo_data_folder="${params.PPI_Coevolution}/CoEvo_data_STRING11.5/"

// move this two parameter here as it could be different for each species, phylum
// can also make them as value channle ?
params.currentSpe_TaxID="511145"
params.current_EggNOG_maxLevel="1224"
        
        
params.input_root_folder="${params.CoEvo_data_folder}allPPI_${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_eggNOGfilteredData"

params.DCA_coevolutoin_path="${params.input_root_folder}/coevolutoin_result_DCA/"


include {prepareSingleMSA_workflow} from "./modules/prepareSingleMSA_workflow.nf"


include {preparePairedMSA_oneRunWithNf90;preparePairedMSA_removeHomologousPairs} from  "./modules/preparePairedMSA.nf"

// include {coevolutionComputation_mfDCA_createdFolder;coevolutionComputation_mfDCA} from "./modules/coevolutionComputation.nf"


workflow {    
    

    println "scriptFile: " + workflow.scriptFile
    println "projectDir: " + workflow.projectDir
    println "launchDir: " + workflow.launchDir
    println "workDir: " + workflow.workDir
    println "configFiles: " + workflow.configFiles
    
    
    prepareSingleMSA_workflow()
        
    // ************Prepare paired MSA **********
    // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
    // here for the because the output generate millions of output,  to prevent completely re-run of progrom in case of accidentaly interuption?
    // some input channle should be the final output path, not the the temperaoy woring directory of previous path 
    // but actully if last step finished correctly , its okay, more the problem of last step compute DCA
    preparePairedMSA_oneRunWithNf90_ch=preparePairedMSA_oneRunWithNf90(prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch.currentSpeMiddleDataPath,prepareSingleMSA_singleMSAGapsFiltering_ch.currentSpe_msa_removeGaps_path)
    preparePairedMSA_removeHomologousPairs_ch=preparePairedMSA_removeHomologousPairs(preparePairedMSA_oneRunWithNf90_ch.pairedMSA_Nf90_csv,preparePairedMSA_oneRunWithNf90_ch.pairedMSA_Nf90_folder)
    
    
}


    
    




/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run PairedMSA_preprocessing_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

to run for other species and phylum
nextflow run PairedMSA_preprocessing_workflow.nf --currentSpe_TaxID "1274374" --current_EggNOG_maxLevel "1239" -params-file wc-params.json -c nextflow.config -resume

Subject_tuple=('1224', '511145')
Query_tupleList=[("1239","1274374"),('201174', '105422'), ('976', '411476'),] 
*/


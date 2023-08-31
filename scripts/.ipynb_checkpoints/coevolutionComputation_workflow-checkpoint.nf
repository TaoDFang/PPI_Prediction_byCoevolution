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
params.currentSpe_TaxID="511145"
params.current_EggNOG_maxLevel="1224"
        
        
params.input_root_folder="${params.CoEvo_data_folder}allPPI_${params.currentSpe_TaxID}_EggNOGmaxLevel${params.current_EggNOG_maxLevel}_eggNOGfilteredData"

params.DCA_coevolutoin_path="${params.input_root_folder}/coevolutoin_result_DCA/"


include {downLoadOtherRawFiles; downLoadRawFastaFile; prepareFastaDataBySpecies;moveOnlyBacteriaSepcies } from './modules/RawFastaFilesAndMetaData'

include {prepareSingleMSA_ParseCurSpeFastaByProteins;prepareSingleMSA_RemoveRedundantProteins;prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs;prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas;prepareSingleMSA_SeedAlignment;prepareSingleMSA_SeedAlignment_filtering;prepareSingleMSA_ProfileHMM_ClustoMA;prepareSingleMSA_ProfileHMM_hmmbuild;prepareSingleMSA_HMMAllAlignment;prepareSingleMSA_singleMSAGapsFiltering} from "./modules/prepareSingleMSA"


include {preparePairedMSA_oneRunWithNf90;preparePairedMSA_removeHomologousPairs} from  "./modules/preparePairedMSA"

include {coevolutionComputation_mfDCA_createdFolder;coevolutionComputation_mfDCA} from "./modules/coevolutionComputation"


workflow {    
    

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
    
    
    
    // ************Prepare single MSA **********
    //http://localhost:8206/lab/workspaces/auto-I/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
    //http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.py

    prepareSingleMSA_ParseCurSpeFastaByProteins_ch=prepareSingleMSA_ParseCurSpeFastaByProteins(prepareFastaDataBySpecies_ch.STRING_fastaBySpecies_Folder)
    
    prepareSingleMSA_RemoveRedundantProteins_ch=prepareSingleMSA_RemoveRedundantProteins(prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpe_fastaData)
    
    
    prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs_ch=prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs(prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpe_fastaData,downLoadOtherRawFiles_ch.eggNOG_folder,downLoadOtherRawFiles_ch.species_file,downLoadOtherRawFiles_ch.species_tree_file)
    
    prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch=prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas(prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs_ch.currentSpe_currentMaxLevel_orthologs,prepareSingleMSA_RemoveRedundantProteins_ch.redundant_proteins_csvFile,moveOnlyBacteriaSepcies_ch.STRING_fastaByBacteriaSpecies_Folder,prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpe_fastaData)
    
    prepareSingleMSA_SeedAlignment_ch=prepareSingleMSA_SeedAlignment(prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpeProSeqPath_ByProteins,prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch.currentSpe_OrthologousGroup_Fa_path)
    
    prepareSingleMSA_SeedAlignment_filtering_ch=prepareSingleMSA_SeedAlignment_filtering(prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpeProSeqPath_ByProteins,moveOnlyBacteriaSepcies_ch.STRING_fastaByBacteriaSpecies_Folder,prepareSingleMSA_SeedAlignment_ch.currentSpe_phmmer_outPath)
    prepareSingleMSA_ProfileHMM_ClustoMA_ch=prepareSingleMSA_ProfileHMM_ClustoMA(prepareSingleMSA_SeedAlignment_filtering_ch.currentSpe_phmmer_OrthologousGroup_path)
    prepareSingleMSA_ProfileHMM_hmmbuild_ch=prepareSingleMSA_ProfileHMM_hmmbuild(prepareSingleMSA_ProfileHMM_ClustoMA_ch.currentSpe_ClustoMSA_path)
    prepareSingleMSA_HMMAllAlignment_ch=prepareSingleMSA_HMMAllAlignment(prepareSingleMSA_ProfileHMM_hmmbuild_ch.currentSpe_hmm_profiles_path,prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch.currentSpe_OrthologousGroup_Fa_path)
    prepareSingleMSA_singleMSAGapsFiltering_ch=prepareSingleMSA_singleMSAGapsFiltering(prepareSingleMSA_HMMAllAlignment_ch.currentSpe_hmm_align_path,prepareSingleMSA_ProfileHMM_ClustoMA_ch.currentSpe_ClustoMSA_path,prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch.currentSpeMiddleDataPath)
    
    
    // ************Prepare paired MSA **********
    // http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
    // here for the because the output generate millions of output,  to prevent completely re-run of progrom in case of accidentaly interuption?
    // some input channle should be the final output path, not the the temperaoy woring directory of previous path 
    // but actully if last step finished correctly , its okay, more the problem of last step compute DCA
    preparePairedMSA_oneRunWithNf90_ch=preparePairedMSA_oneRunWithNf90(prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch.currentSpeMiddleDataPath,prepareSingleMSA_singleMSAGapsFiltering_ch.currentSpe_msa_removeGaps_path)
    preparePairedMSA_removeHomologousPairs_ch=preparePairedMSA_removeHomologousPairs(preparePairedMSA_oneRunWithNf90_ch.pairedMSA_Nf90_csv,preparePairedMSA_oneRunWithNf90_ch.pairedMSA_Nf90_folder)
    
    
    
    // ************Compute DCA  **********
    // !!!!
    //we need to do DCA compuation for millions of PPs , if interupt, better not re-start from scratch 
    // in this case, better not use output to the temperory working directory and then copy to publishDir
    // but direct use final output folder as the input channel or as input parameters 
     coevolutionComputation_mfDCA_createdFolder()
    DCA_coevolutoin_path_ch=Channel.fromPath(params.DCA_coevolutoin_path,type:'dir')
    //below input channle name can put chagned
        coevolutionComputation_mfDCA_ch=coevolutionComputation_mfDCA(DCA_coevolutoin_path_ch,prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch.currentSpeMiddleDataPath,preparePairedMSA_removeHomologousPairs_ch.PPIInfoBeforeCoEvoComp_csv,preparePairedMSA_oneRunWithNf90_ch.pairedMSA_Nf90_folder)
    
    
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
nextflow run PairedMSA_preprocessing_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

*/

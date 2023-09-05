#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// or this one beter ? http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_allPPI_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels.py
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1



//parameter has to be define before include state to be able be used by the external modules 
params.newSTRING_rootFolder="${params.PPI_Coevolution}/STRING_data_11.5" //this folder is just newSTRING_rootFolder


// here need to define  species and phylum params as valuen channel
// as it could be different to the same prepareSingleMSA_workflow but with different species and phylum params as input 




include {RawFastaFilesAndMetaData_workflow } from './RawFastaFilesAndMetaData_workflow.nf'

include {prepareSingleMSA_ParseCurSpeFastaByProteins;prepareSingleMSA_RemoveRedundantProteins;prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs;prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas;prepareSingleMSA_SeedAlignment;prepareSingleMSA_SeedAlignment_filtering;prepareSingleMSA_ProfileHMM_ClustoMA;prepareSingleMSA_ProfileHMM_hmmbuild;prepareSingleMSA_HMMAllAlignment;prepareSingleMSA_singleMSAGapsFiltering} from "./modules/prepareSingleMSA"




workflow prepareSingleMSA_workflow{    
    
    take: 
        
        currentSpe_TaxID_ch
        current_EggNOG_maxLevel_ch
        STRING_fastaBySpecies_Folder
        eggNOG_folder
        species_file
        species_tree_file
        STRING_fastaByBacteriaSpecies_Folder
    main:

        println "scriptFile: " + workflow.scriptFile
        println "projectDir: " + workflow.projectDir
        println "launchDir: " + workflow.launchDir
        println "workDir: " + workflow.workDir
        println "configFiles: " + workflow.configFiles
    
        // ************Prepare single MSA **********
        //http://localhost:8206/lab/workspaces/auto-I/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
        //http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/script_CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.py

        prepareSingleMSA_ParseCurSpeFastaByProteins_ch=prepareSingleMSA_ParseCurSpeFastaByProteins(currentSpe_TaxID_ch,STRING_fastaBySpecies_Folder)

        prepareSingleMSA_RemoveRedundantProteins_ch=prepareSingleMSA_RemoveRedundantProteins(currentSpe_TaxID_ch,prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpe_fastaData)


        prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs_ch=prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs(currentSpe_TaxID_ch,current_EggNOG_maxLevel_ch,prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpe_fastaData,eggNOG_folder,species_file,species_tree_file)

        prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch=prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas(currentSpe_TaxID_ch,current_EggNOG_maxLevel_ch,prepareSingleMSA_PreprocessEggnogOrthologGroup_chooseOrthologs_ch.currentSpe_currentMaxLevel_orthologs,prepareSingleMSA_RemoveRedundantProteins_ch.redundant_proteins_csvFile,STRING_fastaByBacteriaSpecies_Folder,prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpe_fastaData)

        prepareSingleMSA_SeedAlignment_ch=prepareSingleMSA_SeedAlignment(currentSpe_TaxID_ch,current_EggNOG_maxLevel_ch,prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpeProSeqPath_ByProteins,prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch.currentSpe_OrthologousGroup_Fa_path)

        prepareSingleMSA_SeedAlignment_filtering_ch=prepareSingleMSA_SeedAlignment_filtering(currentSpe_TaxID_ch,current_EggNOG_maxLevel_ch,prepareSingleMSA_ParseCurSpeFastaByProteins_ch.currentSpeProSeqPath_ByProteins,STRING_fastaByBacteriaSpecies_Folder,prepareSingleMSA_SeedAlignment_ch.currentSpe_phmmer_outPath)
        prepareSingleMSA_ProfileHMM_ClustoMA_ch=prepareSingleMSA_ProfileHMM_ClustoMA(currentSpe_TaxID_ch,current_EggNOG_maxLevel_ch,prepareSingleMSA_SeedAlignment_filtering_ch.currentSpe_phmmer_OrthologousGroup_path)
        prepareSingleMSA_ProfileHMM_hmmbuild_ch=prepareSingleMSA_ProfileHMM_hmmbuild(currentSpe_TaxID_ch,current_EggNOG_maxLevel_ch,prepareSingleMSA_ProfileHMM_ClustoMA_ch.currentSpe_ClustoMSA_path)
        prepareSingleMSA_HMMAllAlignment_ch=prepareSingleMSA_HMMAllAlignment(currentSpe_TaxID_ch,current_EggNOG_maxLevel_ch,prepareSingleMSA_ProfileHMM_hmmbuild_ch.currentSpe_hmm_profiles_path,prepareSingleMSA_PreprocessEggnogOrthologGroup_collectingOGFastas_ch.currentSpe_OrthologousGroup_Fa_path)
        prepareSingleMSA_singleMSAGapsFiltering_ch=prepareSingleMSA_singleMSAGapsFiltering(currentSpe_TaxID_ch,current_EggNOG_maxLevel_ch,prepareSingleMSA_HMMAllAlignment_ch.currentSpe_hmm_align_path,prepareSingleMSA_ProfileHMM_ClustoMA_ch.currentSpe_ClustoMSA_path)
    
    emit:
        currentSpeMSAGapsFilteringMetaFolder=prepareSingleMSA_singleMSAGapsFiltering_ch.currentSpeMSAGapsFilteringMetaFolder
        currentSpe_msa_removeGaps_path=prepareSingleMSA_singleMSAGapsFiltering_ch.currentSpe_msa_removeGaps_path
        
    
    
}


workflow {
        // ************Download all metadata and sequencing data**********
        // here better to seprate subworkflow RawFastaFilesAndMetaData_workflow prepareSingleMSA_workflow
        // as these two workflow are kind of independent and first workflow provide output to the next workflow 
        //also no need to use contional run now, as it will be excuted once here anyway
        //here to check if RawFastaFilesAndMetaData_workflow process has been run before ; in worse case, use a params to decied to re run or not 
        // downloadedTreeFile = file("${params.RawData_Folder}/species.tree.v11.5.txt") // interesting ,here has to be double quoto""
        // // println "${myFile.getParent()};${myFile.getName()};${myFile.exists()}"
        // if (downloadedTreeFile.exists()) {
        //     println "RawFastaFilesAndMetaData_workflow has been run before,skip it "
        // }
        // else {
        //     RawFastaFilesAndMetaData_workflow() 
        // }
    
        RawFastaFilesAndMetaData_workflow() 
    
        prepareSingleMSA_workflow(Channel.value("${params.query_currentSpe_TaxID}"),
                             Channel.value("${params.query_current_EggNOG_maxLevel}"),
                                  RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                                  RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                                  RawFastaFilesAndMetaData_workflow.out.species_file,
                                  RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                                  RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                                 )
    
        prepareSingleMSA_workflow(Channel.value("${params.subject1_currentSpe_TaxID}"),
                             Channel.value("${params.subject1_current_EggNOG_maxLevel}"),
                                  RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                                  RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                                  RawFastaFilesAndMetaData_workflow.out.species_file,
                                  RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                                  RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                                 )
    
        prepareSingleMSA_workflow(Channel.value("${params.subject2_currentSpe_TaxID}"),
                             Channel.value("${params.subject2_current_EggNOG_maxLevel}"),
                                  RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                                  RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                                  RawFastaFilesAndMetaData_workflow.out.species_file,
                                  RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                                  RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                                 )
    
        prepareSingleMSA_workflow(Channel.value("${params.subject3_currentSpe_TaxID}"),
                             Channel.value("${params.subject3_current_EggNOG_maxLevel}"),
                                  RawFastaFilesAndMetaData_workflow.out.STRING_fastaBySpecies_Folder,
                                  RawFastaFilesAndMetaData_workflow.out.eggNOG_folder,
                                  RawFastaFilesAndMetaData_workflow.out.species_file,
                                  RawFastaFilesAndMetaData_workflow.out.species_tree_file,
                                  RawFastaFilesAndMetaData_workflow.out.STRING_fastaByBacteriaSpecies_Folder,
                                 )
    
}
    



/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training

when this folder is in /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run prepareSingleMSA_workflow.nf  -params-file wc-params.json -c nextflow.config -resume


with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

*/


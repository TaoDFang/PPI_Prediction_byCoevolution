#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1



params.RawData_Folder = '/mnt/mnemo6/tao/nextflow/STRING_Data_11.5/'
//STRING_root_folder="/mnt/mnemo6/tao/STRING_Data_11.5/"





workflow {
    Test()
    Test.out.view()
    
    //somehow this chnnel didnt give me folder path, but keep give me "Base_Folder: unbound variable" error, maybe its better for the files ?, and if remove "type: 'dir'", do nothing 
    // Problem solved by replace ''' ''' to """ """" for "script" section in downLoadRawFastaFile process
    // RawData_Folder_ch=Channel.fromPath('${params.RawData_Folder}/*',type:'dir')
    RawData_Folder_ch=Channel.fromPath(params.RawData_Folder,type:'dir')
    // here cause error : Aug-14 16:32:42.550 [Actor Thread 12] DEBUG nextflow.util.CacheHelper - Unable to get file attributes file: /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts/${params.RawData_Folder}/* -- Cause: java.nio.file.NoSuchFileException: /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts/${params.RawData_Folder}/*
    // reason: usuage of ${parameter} is for Script parameters, not use in nextflow other part 
    RawData_Folder_ch.view()
    
    downLoadRawFastaFile_ch=downLoadRawFastaFile(RawData_Folder_ch)
    //downLoadRawFastaFile.out.view()
    downLoadRawFastaFile_ch.view()
}

    
process Test {
    
    label "simple_process"
    
    output:
    stdout
        
    script:
    '''
    echo test started
    echo test finished
    '''
}


//http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/PrepareFastaDataBySpecies.ipynb
// download all fasta seq data from new string website 
process downLoadRawFastaFile {
    
    label "simple_process"
    
    input: 
    path i_RawData_Folder
    
    output:
    //stdout
    // path "${i_RawData_Folder}/STRINGSequencesBySpecies/", emit: STRING_fastaBySpecies_Folder
    // path "${i_RawData_Folder}/STRINGBacteriaSequencesBySpecies/", emit: STRING_fastaByBacteriaSpecies_Folder
    val "${STRING_fastaBySpecies_Folder}"
    val "${STRING_fastaByBacteriaSpecies_Folder}"
        
    
    script:
    """
    
    echo process downLoadRawFastaFile started
    
    echo ${i_RawData_Folder} #output: STRING_Data_11.5
    echo ${params.RawData_Folder} #output: /mnt/mnemo6/tao/nextflow/STRING_Data_11.5/
    #my understanding is the channel is more  important for the "input" of intermediate task, not the first task 
    
    
    STRING_fastaBySpecies_Folder="${i_RawData_Folder}/STRINGSequencesBySpecies/"
    STRING_fastaByBacteriaSpecies_Folder="${i_RawData_Folder}/STRINGBacteriaSequencesBySpecies/"
    
    #here notice for bash variables(defined withing script block), it need to be used by "\" charater
    echo \${STRING_fastaBySpecies_Folder}
    echo \${STRING_fastaByBacteriaSpecies_Folder}
    
    

    
    echo process downLoadRawFastaFile finished
    
    """
}



//     # where here generated multiple layer folder ??
//     cd ${i_RawData_Folder}
//     wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz -P .
//     gzip -d protein.sequences.v11.5.fa.gz
    


// process prepareFastaDataBySpecies {
//     here run a python file, input is output of last process 
//     script:
//     """
//     STRING_fastaBySpecies_Folder="${i_RawData_Folder}/STRINGSequencesBySpecies/"
//     STRING_fastaByBacteriaSpecies_Folder="${i_RawData_Folder}/STRINGBacteriaSequencesBySpecies/"
    
//     #here notice for bash variables(defined withing script block), it need to be used by "\" charater
//     echo \${STRING_fastaBySpecies_Folder}
//     echo \${STRING_fastaByBacteriaSpecies_Folder}


    
//     """
// }

/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run PairedMSA_preprocessing_workflow.nf -params-file wc-params.json -c nextflow.config -resume
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

    # download all fasta seq data from new string website 
    mkdir -p ${Base_Folder}/STRING_Data_11.5
    cd ${Base_Folder}/STRING_Data_11.5
    wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz .
    gzip -d protein.sequences.v11.5.fa.gz
    

*/

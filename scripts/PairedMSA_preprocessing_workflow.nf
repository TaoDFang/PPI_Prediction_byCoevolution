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
    //somehow this chnnel didnt give me folder path, but keep give me "Base_Folder: unbound variable" error, maybe its better for the files ?, and if remove "type: 'dir'", do nothing 
    // Problem solved by replace ''' ''' to """ """" for "script" section in downLoadRawFastaFile process
    // RawData_Folder_ch=Channel.fromPath('${params.RawData_Folder}/*',type:'dir')
    RawData_Folder_ch=Channel.fromPath(params.RawData_Folder,type:'dir')
    // here cause error : Aug-14 16:32:42.550 [Actor Thread 12] DEBUG nextflow.util.CacheHelper - Unable to get file attributes file: /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts/${params.RawData_Folder}/* -- Cause: java.nio.file.NoSuchFileException: /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts/${params.RawData_Folder}/*
    // reason: usuage of ${parameter} is for Script parameters, not use in nextflow other part 
    RawData_Folder_ch.view()
    
    downLoadRawFastaFile_ch=downLoadRawFastaFile(RawData_Folder_ch)
    //downLoadRawFastaFile.out.view()
    //downLoadRawFastaFile_ch.rawFasta_file.view()
    // downLoadRawFastaFile_ch.test_output.view()
    //downLoadRawFastaFile_ch.view(), this is not working when there are multiple output files that are not from same object, alternatively use tuple
}



    

//http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/PrepareFastaDataBySpecies.ipynb
// download all fasta seq data from new string website 
process downLoadRawFastaFile {
    
    label "simple_process"
    
    input: 
    path i_RawData_Folder
    
    
    //when: somehow this process is repeatly executed with -resume option, i think the resason is that  folder "i_RawData_Folder" changed with the output?
    // the problem is of using when is that the output channel will be empty, this will cause problem from next process ß
        
    output:
    //stdout
        
    //here seem path has to be the output from some where in the script , 
    //path "511145.protein.sequences.v11.5.fa", emit: rawFasta_file
    #path "process_finished.txt", emit: test_output

    
    script:
    """
    
    echo process downLoadRawFastaFile started 
    
    echo ${i_RawData_Folder} #output: STRING_Data_11.5
    echo ${params.RawData_Folder}  #output: /mnt/mnemo6/tao/nextflow/STRING_Data_11.5/
    #my understanding is the channel is more  important for the "input" of intermediate task, not the first task 


    #cd ${i_RawData_Folder} #here this command cause error, no output, why ???, it seem create a new folder inside folder i_RawData_Folder ?
    #wget https://stringdb-downloads.org/download/protein.sequences.v11.5/511145.protein.sequences.v11.5.fa.gz -P ${i_RawData_Folder}
    
    #gunzip -c "${i_RawData_Folder}/511145.protein.sequences.v11.5.fa.gz" > 511145.protein.sequences.v11.5.fa
    
    echo process downLoadRawFastaFile finished  # this one is cached , 
    #echo process downLoadRawFastaFile finished  >process_finished.txt # this one is not cached, why ??  because it also change timestamp of folder "i_RawData_Folder"?

    
    """
}




//     # why here generated multiple layer folder ??
//     cd ${i_RawData_Folder}
//     wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz -P .
//     gzip -d protein.sequences.v11.5.fa.gz > protein.sequences.v11.5.fa
    


//     STRING_fastaBySpecies_Folder="${i_RawData_Folder}/STRINGSequencesBySpecies/"
//     STRING_fastaByBacteriaSpecies_Folder="${i_RawData_Folder}/STRINGBacteriaSequencesBySpecies/"
//     #here notice for bash variables(defined withing script block), it need to be used by "\" charater
//     echo \${STRING_fastaBySpecies_Folder} >
//     echo \${STRING_fastaByBacteriaSpecies_Folder}
    
    


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
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

    # download all fasta seq data from new string website 
    mkdir -p ${Base_Folder}/STRING_Data_11.5
    cd ${Base_Folder}/STRING_Data_11.5
    wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz .
    gzip -d protein.sequences.v11.5.fa.gz
    

*/

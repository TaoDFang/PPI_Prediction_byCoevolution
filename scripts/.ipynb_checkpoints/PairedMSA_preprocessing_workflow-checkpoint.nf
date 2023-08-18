#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 

// for the text editor format , choose groovy 

//https://github.com/TaoDFang/MNF/issues/98
//https://github.com/TaoDFang/PPI_Prediction_byCoevolution/issues/1



params.RawData_Folder = '/mnt/mnemo6/tao/nextflow/STRING_Data_11.5'
params.STRING_fastaBySpecies_Folder="${params.RawData_Folder}/STRINGSequencesBySpecies/"
//STRING_root_folder="/mnt/mnemo6/tao/STRING_Data_11.5/"





workflow {    
    //somehow this chnnel didnt give me folder path, but keep give me "Base_Folder: unbound variable" error, maybe its better for the files ?, and if remove "type: 'dir'", do nothing 
    // Problem solved by replace ''' ''' to """ """" for "script" section in downLoadRawFastaFile process
    // RawData_Folder_ch=Channel.fromPath('${params.RawData_Folder}/*',type:'dir')
    //RawData_Folder_ch=Channel.fromPath(params.RawData_Folder,type:'dir')
    // here cause error : Aug-14 16:32:42.550 [Actor Thread 12] DEBUG nextflow.util.CacheHelper - Unable to get file attributes file: /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts/${params.RawData_Folder}/* -- Cause: java.nio.file.NoSuchFileException: /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts/${params.RawData_Folder}/*
    // reason: usuage of ${parameter} is for Script parameters, not use in nextflow other part 
    //RawData_Folder_ch.view()
    
    //downLoadRawFastaFile_ch=downLoadRawFastaFile(RawData_Folder_ch)
    downLoadRawFastaFile_ch=downLoadRawFastaFile()
    //downLoadRawFastaFile.out.view()
    // downLoadRawFastaFile_ch.rawFasta_file.view()
    // downLoadRawFastaFile_ch.test_output.view()
    //downLoadRawFastaFile_ch.view(), this is not working when there are multiple output files that are not from same object, alternatively use tuple
    
    
    prepareFastaDataBySpecies_ch=prepareFastaDataBySpecies(downLoadRawFastaFile_ch.rawFasta_file)
    //prepareFastaDataBySpecies.out.view()
    
    // here use channel not params.STRING_fastaBySpecies_Folder directly to triger next process ""
    //STRING_fastaBySpecies_Folder_ch=Channel.fromPath(params.STRING_fastaBySpecies_Folder,type:'dir')
    moveOnlyBacteriaSepcies(prepareFastaDataBySpecies_ch)
    moveOnlyBacteriaSepcies.out.view()
    
    
}



    

//http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/PrepareFastaDataBySpecies.ipynb
// download all fasta seq data from new string website 
process downLoadRawFastaFile {
    
    label "simple_process"
    
    
    publishDir "${params.outdir}"
    
    
    // input: 
    // path i_RawData_Folder
    
    
    //when: somehow this process is repeatly executed with -resume option, i think the resason is that  folder "i_RawData_Folder" changed with the output?
    // the problem is of using when is that the output channel will be empty, this will cause problem from next process ß
        
    output:
    //stdout
        
    //here seem path has to be the output from somewhere in the script , 
    path "protein.sequences.v11.5.fa", emit: rawFasta_file
    path "species.v11.5.txt", emit: STR_species_mem
    path "process_finished.txt", emit: test_output

    
    script:
    """
    
    echo process downLoadRawFastaFile started 
    
    echo ${params.RawData_Folder}  # output: /mnt/mnemo6/tao/nextflow/STRING_Data_11.5/

    # cd \${RawData_Folder} here this command cause error, no output, why ???, it seem create a new folder inside folder RawData_Folder ?
    
    wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz  -P ${params.RawData_Folder}
    gunzip -c "${params.RawData_Folder}/protein.sequences.v11.5.fa.gz" >  protein.sequences.v11.5.fa
    
    # also download other STRING meatafile 
    wget https://stringdb-downloads.org/download/species.v11.5.txt -P ${params.RawData_Folder} -O species.v11.5.txt # here has to use -O to get actually output

    echo process downLoadRawFastaFile finished  > process_finished.txt 
    
    """
}




    

//     STRING_fastaBySpecies_Folder="${i_RawData_Folder}/STRINGSequencesBySpecies/"
//     STRING_fastaByBacteriaSpecies_Folder="${i_RawData_Folder}/STRINGBacteriaSequencesBySpecies/"
//     #here notice for bash variables(defined withing script block), it need to be used by "\" charater
//     echo \${STRING_fastaBySpecies_Folder} >
//     echo \${STRING_fastaByBacteriaSpecies_Folder}
    
    


process prepareFastaDataBySpecies {
    tag "process  prepareFastaDataBySpecies"
    publishDir "${params.outdir}"
    // cpus 32
    // memory 100.GB
    
    conda "/mnt/mnemo5/tao/anaconda3/envs/ipykernel_py3"  // specify conda enviroment here works, but why in configuration file not working 
    debug true //echo true echo directive is deprecated 
    
    input:
        path rawFasta_file
        //tuple(path(rawFasta_file),path(STR_species_mem))

    
    output: 
       path "STRINGSequencesBySpecies", type : "dir", emit: STRING_fastaBySpecies_Folder
       //here could not use full path,  "${params.RawData_Folder}/STRINGSequencesBySpecies/",
        //other wise get error : file `/mnt/mnemo6/tao/nextflow/STRING_Data_11.5/STRINGSequencesBySpecies/` is outside the scope of the process work directory: /mnt/mnemo6/tao/nextflow/work/fd/ceda6c2edc10b2d7971b44198ce727
    
    
        // env to cappute variable defined in script  ,but this step still created intemediated files 
        //https://stackoverflow.com/questions/70862345/how-to-declare-a-new-nextflow-variable-in-a-shell-block
        //env STRING_fastaBySpecies_Folder, but seem it can out used as path channel later, but literally a string "STRING_fastaBySpecies_Folder"
        
    //     path "${params.RawData_Folder}/STRINGSequencesBySpecies/", emit: STRING_fastaBySpecies_Folder, this one also not working  due the to ${params.RawData_Folder}/STRINGSequencesBySpecies/ is not the in thge working directory of this folder 
     //   stdout
    //     path "STRING_fastaBySpecies_Folder_path", emit: STRING_fastaBySpecies_Folder
    // if the output is a folder to be used later, can used channel defined outside the process ?
    // or could use stdout
    
    script:
        
    """
    echo ${rawFasta_file}
    echo $CONDA_DEFAULT_ENV  # this is not specified conda enviroment, but seem bewlowing python use the correct one 
    

    
    #here nto use params options diretly to be able to  trigle next process  moveOnlyBacteriaSepcies
    STRING_fastaBySpecies_Folder_tmp="${params.RawData_Folder}/STRINGSequencesBySpecies_tmp/"
    mkdir -p  \${STRING_fastaBySpecies_Folder_tmp}
    
    # here use without python instead of 'python prepareFastaDataBySpecies.py', its nextflow special : https://carpentries-incubator.github.io/workflows-nextflow/aio/index.html
    # here in conda enviroment "nf-training", python3 works but not python
    python ${projectDir}/python_scripts/prepareFastaDataBySpecies.py --rawFasta_file ${rawFasta_file} --STRING_fastaBySpecie \${STRING_fastaBySpecies_Folder_tmp}
    
    mv \${STRING_fastaBySpecies_Folder_tmp} "STRINGSequencesBySpecies" # use this to be able to export folder as process output 
    # altinatively, use solution output folder name  (better form python ) to a file and read content from the file 
    
    """
}


process moveOnlyBacteriaSepcies {
    input: 
        path x
    
    output:
        stdout
    
    script:
    """
        echo process moveOnlyBacteriaSepcies
        echo ${x}

    """
    
}

// the thing is next copy process should depend on successfuly exeution of last command ?, can output something folder from last script ??nex
// still need to use STRING_fastaBySpecies_Folder as channel 


//     #here notice for bash variables(defined withing script block), it need to be used by "\" charater
//     #STRING_fastaBySpecies_Folder="${params.RawData_Folder}/STRINGSequencesBySpecies/"
//     #echo \${STRING_fastaBySpecies_Folder}
//     #echo \${STRING_fastaBySpecies_Folder} > "STRING_fastaBySpecies_Folder_path", not actully working to output folder STRING_fastaBySpecies_Folder
    




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

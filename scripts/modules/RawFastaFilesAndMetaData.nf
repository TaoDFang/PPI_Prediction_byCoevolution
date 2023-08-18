// ************Download all metadata and sequencing data**********

//http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/PrepareFastaDataBySpecies.ipynb
// download all fasta seq data from new string website 
process downLoadRawFastaFile {
    
    label "simple_process"
    
    
    publishDir "${params.outdir}"
    
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
    
    # when using full path, the output is not in process working directory and thus no need to use publishDir directive
    wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz  -P ${params.RawData_Folder}
    gunzip -c "${params.RawData_Folder}/protein.sequences.v11.5.fa.gz" >  protein.sequences.v11.5.fa
    
    # also download other STRING meatafile 
    wget https://stringdb-downloads.org/download/species.v11.5.txt -P ${params.RawData_Folder} -O species.v11.5.txt # here has to use -O to get actually output

    echo process downLoadRawFastaFile finished  > process_finished.txt 
    
    """
}


process prepareFastaDataBySpecies {
    tag "process  prepareFastaDataBySpecies"
    publishDir "${params.outdir}" , mode: "copy"
    // cpus 32
    // memory 100.GB
    
    conda "/mnt/mnemo5/tao/anaconda3/envs/ipykernel_py3"  // specify conda enviroment here works, but why in configuration file not working , works now 
    debug true //echo true echo directive is deprecated 
    
    input:
        path rawFasta_file
        //tuple(path(rawFasta_file),path(STR_species_mem))
    
    output: 
       path "STRINGSequencesBySpecies/", type : "dir", emit: STRING_fastaBySpecies_Folder

    script:
        
    """
    echo ${rawFasta_file}
    echo $CONDA_DEFAULT_ENV  # this is not specified  correct conda enviroment, but seem bewlowing python use the correct one 
    
    #here do not  use params options diretly to be able to  trigle next process  moveOnlyBacteriaSepcies
    STRING_fastaBySpecies_Folder_tmp="${params.RawData_Folder}/STRINGSequencesBySpecies_tmp/"
    mkdir -p  \${STRING_fastaBySpecies_Folder_tmp}
    
    python ${projectDir}/python_scripts/prepareFastaDataBySpecies.py --rawFasta_file ${rawFasta_file} --STRING_fastaBySpecie \${STRING_fastaBySpecies_Folder_tmp}
    
    mv \${STRING_fastaBySpecies_Folder_tmp} "STRINGSequencesBySpecies/" # use this to be able to export folder as process output 
    # altinatively, use solution output folder name  (better form python ) to a file and read content from the file 
    
    """
}


process moveOnlyBacteriaSepcies {
    publishDir "${params.outdir}", mode: "copy"
    
    input: 
        path STR_species_mem
        path STRING_fastaBySpecie

    
    output:
        path "STRINGBacteriaSequencesBySpecies/", type: "dir", emit: STRING_fastaByBacteriaSpecies_Folder
    
    script:
    """
        echo process moveOnlyBacteriaSepcies
        mkdir -p "STRINGBacteriaSequencesBySpecies/"
        
        # when using relative path, the output is not in working directory and thus no need need to use  publishDir directive to actulyl copy the data from/out of process working directory 
        python ${projectDir}/python_scripts/moveOnlyBacteriaSepcies.py --STR_species_mem_file ${STR_species_mem} --STRING_fastaBySpecie ${STRING_fastaBySpecie} --STRING_fastaByBacteriaSpecies "STRINGBacteriaSequencesBySpecies/"
        
    """
    
}

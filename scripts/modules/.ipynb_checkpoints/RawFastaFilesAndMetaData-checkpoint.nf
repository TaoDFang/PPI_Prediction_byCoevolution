// ************Download all metadata and sequencing data**********



process test_configuration {
    
    label "simple_process"
    
    
    publishDir "${params.RawData_Folder}", mode: "copy" 
    debug true //echo true echo directive is depreca
    
    output:
        path "process_finished.txt", emit: test_output

    script:
    """
        echo ${params.RawData_Folder} >> process_finished.txt 
        echo $PATH >> process_finished.txt 
        echo $CONDA_DEFAULT_ENV >> process_finished.txt  # here problem is to print conda envs of nf-training, but of the current process
        conda info --envs >> process_finished.txt 
        python --version >> process_finished.txt  # for nf-training enviroment, python doest work
        echo ${params.container_path} == '' >> process_finished.txt
        echo ${workflow.containerEngine } >> process_finished.txt
        echo ${workflow.container }  >> process_finished.txt
        hmmbuild -h  >> process_finished.txt  # this works only for conda envs sequence_tools_conda
        echo process test_configuration finished  >> process_finished.txt 
        echo process test_configuration finished  >> process_finished.txt 

    """
}

process downLoadOtherRawFiles {
    
    label "simple_process"
    
    
    publishDir "${params.RawData_Folder}", mode: "copy" 
    debug true //echo true echo directive is depreca
    
    output:
    //stdout

    path "species.v11.5.txt", type: "file", emit: species_file
    path "species.tree.v11.5.txt", type: "file", emit: species_tree_file
    path "eggnog5AddSTRING11.5_Species/", type: "dir", emit: eggNOG_folder 

    
    script:
    """
        echo ${params.RawData_Folder}
        echo $CONDA_DEFAULT_ENV
        
        wget https://stringdb-downloads.org/download/species.v11.5.txt -O species.v11.5.txt  # here -N is not necesseary since its in temperaroy processing working directoy ,  -N
        wget https://stringdb-downloads.org/download/species.tree.v11.5.txt  -O species.tree.v11.5.txt 
        
        #download eggnog file
        eggNOG_folder="eggnog5AddSTRING11.5_Species/"
        mkdir -p \${eggNOG_folder} # here -p is not necesseary since its in temperaroy processing working directoy 
        wget https://zenodo.org/record/8279323/files/eggnog5AddSTRING11.5_Species.tar.gz?download=1  -O eggnog5AddSTRING11.5_Species.tar.gz
        tar -xf "eggnog5AddSTRING11.5_Species.tar.gz" -C \${eggNOG_folder} #https://linuxhint.com/solve-gzip-stdin-not-gzip-format-error/. without -v to avoid outut infor -v
    """
}







process downLoadRawFastaFile {
    
    label "simple_process"
    
    
    publishDir "${params.RawData_Folder}", mode: "copy"
    
    output:
    //stdout
        
    //here seem path has to be the output from somewhere in the script , 
    path "protein.sequences.v11.5.fa", emit: rawFasta_file


    
    script:
    """
    
    echo process downLoadRawFastaFile started 
    echo  $CONDA_DEFAULT_ENV
    
    echo ${params.RawData_Folder}  # output: /mnt/mnemo6/tao/nextflow/STRING_Data_11.5/

    # cd \${RawData_Folder} here this command cause error, no output, why ???, it seem create a new folder inside folder RawData_Folder ?
    
    
    wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz -O protein.sequences.v11.5.fa.gz  # here actull no need to use  -P \${params.RawData_Folder}, since we will copy data from current process working directory to the output channel
    # when using full path like -P \${params.RawData_Folder}, the output is not in process working directory and thus no need to use publishDir directive, but better not do it 
    gunzip -c protein.sequences.v11.5.fa.gz >  protein.sequences.v11.5.fa
    

    
    """
}


process prepareFastaDataBySpecies {
    tag "process  prepareFastaDataBySpecies"
    publishDir "${params.RawData_Folder}" , mode: "copy"
    
    label "simple_process"
    debug true //echo true echo directive is deprecated 
    
    input:
        path rawFasta_file
        //tuple(path(rawFasta_file),path(STR_species_mem))
    
    output: 
       path "STRINGSequencesBySpecies/", type : "dir", emit: STRING_fastaBySpecies_Folder

    script:
        
    """
    echo  ${rawFasta_file}
    echo  $CONDA_DEFAULT_ENV  # this is not specified  correct conda enviroment, but seem bewlowing python use the correct one 
    
    #here do not  use params options diretly to be able to  trigle next process  moveOnlyBacteriaSepcies
    STRING_fastaBySpecies_Folder="STRINGSequencesBySpecies/"
    mkdir -p  \${STRING_fastaBySpecies_Folder}  # here -p is not necesseary since its in temperaroy processing working directoy 
    python ${projectDir}/python_scripts/prepareFastaDataBySpecies.py --rawFasta_file ${rawFasta_file} --STRING_fastaBySpecie \${STRING_fastaBySpecies_Folder}
    
    """
}


process moveOnlyBacteriaSepcies {
    publishDir "${params.RawData_Folder}", mode: "copy"
    label "simple_process"
    
    
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

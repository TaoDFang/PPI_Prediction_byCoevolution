nextflow.enable.dsl=2

process PrepareFastaDataBySpecies {
    input: 
    path Base_Folder
    
    output:
    stdout
        
    script:
    '''
    
     echo started
    
    '''
}

workflow {
    Base_Folder_ch=Channel.fromPath("/mnt/mnemo6/tao/nextflow/*")
    PrepareFastaDataBySpecies(Base_Folder_ch)
    PrepareFastaDataBySpecies.out.view()
}




//process_publishDir.nf
nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 



// params.Base_Folder = "/mnt/mnemo6/tao/nextflow"
params.Base_Folder = "/mnt/mnemo6/tao/nextflow/text.text"

//http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/PrepareFastaDataBySpecies.ipynb



workflow {
    Test()
    Test.out.view()
    
    //somehow this chnnel didnt give me folder path, but keep give me "Base_Folder: unbound variable" error, maybe its better for the files ?, and if remove "type: 'dir'", do nothing 
    // Base_Folder_ch=Channel.fromPath("${params.Base_Folder}/text.text")
    Base_Folder_ch=Channel.fromPath(params.Base_Folder)
    PrepareFastaDataBySpecies(Base_Folder_ch)
    
    //PrepareFastaDataBySpecies(params.Base_Folder)
    PrepareFastaDataBySpecies.out.view()
}


process PrepareFastaDataBySpecies {
    input: 
    path i_Base_Folder
    
    output:
    stdout
        
    script:
    '''
    echo started
    echo finished
    echo ${i_Base_Folder}
    '''
}


process Test {
    
    output:
    stdout
        
    script:
    '''
    echo test started
    echo test finished
    '''
}


/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run PairedMSA_preprocessing_workflow.nf

    # download all fasta seq data from new string website 
    mkdir -p ${Base_Folder}/STRING_Data_11.5
    cd ${Base_Folder}/STRING_Data_11.5
    wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz .
    gzip -d protein.sequences.v11.5.fa.gz
    

*/

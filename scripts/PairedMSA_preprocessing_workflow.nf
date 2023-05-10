//process_publishDir.nf
nextflow.enable.dsl=2

// the main logic is from http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/CoEvo_EggNOG_preprocessing_STRING1105_varyEggNOGMaxLevels_prepareSTRINPhyPPIBenchmark.ipynb
// compare with one jupyter notebook , the snakemake has better pipeline, (running )enviroment managements 



params.Base_Folder = "/mnt/mnemo6/tao/nextflow/"

//http://localhost:8206/lab/workspaces/auto-Z/tree/code/MNF/notebooks/STRING_Data_11.5/PrepareFastaDataBySpecies.ipynb
process PrepareFastaDataBySpecies {
    input: 
    path Base_Folder
        
    script:
    '''
    mkdir -p ${Base_Folder}nextflow
    
    # download all fasta seq data from new string website 
    mkdir -p ${Base_Folder}nextflow/STRING_Data_11.5
    cd ${Base_Folder}nextflow/STRING_Data_11.5
    wget https://stringdb-static.org/download/protein.sequences.v11.5.fa.gz .
    gzip -d protein.sequences.v11.5.fa.gz
    '''
}

workflow {
    Base_Folder_ch=Channel.fromPath(params.Base_Folder)
    PrepareFastaDataBySpecies(Base_Folder_ch)
}


/*
* optional: test in a tmux sesssion:  tmux attach -t tmux-nextflow 
conda activate nf-training
cd /mnt/mnemo5/tao/PPI_Prediction_byCoevolution/scripts
nextflow run PairedMSA_preprocessing_workflow.nf
*/

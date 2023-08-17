    echo process downLoadRawFastaFile starteFd 
    
    #echo ${i_RawData_Folder} #output: STRING_Data_11.5, use it in this process has a problem, its timestamp changed with the new output 
    echo ${params.RawData_Folder}  #output: /mnt/mnemo6/tao/nextflow/STRING_Data_11.5/
    #my understanding is the channel is more  important for the "input" of intermediate task, not the first task 


    #cd ${i_RawData_Folder} #here this command cause error, no output, why ???, it seem create a new folder inside folder i_RawData_Folder ?
    #wget https://stringdb-downloads.org/download/protein.sequences.v11.5/511145.protein.sequences.v11.5.fa.gz -P ${i_RawData_Folder} #this one is not cached, why ??  because it also change timestamp of folder "i_RawData_Folder"?
    
    wget https://stringdb-downloads.org/download/protein.sequences.v11.5/511145.protein.sequences.v11.5.fa.gz -P ${params.RawData_Folder}}
    
    #gunzip -c "${i_RawData_Folder}/511145.protein.sequences.v11.5.fa.gz" > 511145.protein.sequences.v11.5.fa
    
    echo process downLoadRawFastaFile finished  # this one is cached , 
    echo process downLoadRawFastaFile finished  >process_finished.txt # this one is not cached, why ??  because it also change timestamp of folder "i_RawData_Folder"?

    
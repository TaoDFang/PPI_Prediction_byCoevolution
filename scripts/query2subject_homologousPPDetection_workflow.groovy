#!/usr/bin/env nextflow

nextflow.enable.dsl=2



params.newSTRING_rootFolder="${params.PPI_Coevolution}/STRING_data_11.5" //this folder is just newSTRING_rootFolder
params.homologous_ppPath="${params.PPI_Coevolution}/STRING_data_11.5/homologous_pp" 




// the main logic is from http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_homologousPPDetection.ipynb
// and/or mainly from  http://localhost:8206/lab/workspaces/auto-j/tree/code/MNF/notebooks/STRING_Data_11.5/test_phylumeffect_homologousPPDetection_STRINGPhyBalancePhyla.ipynb


workflow query2subject_homologousPPDetection_workflow{    
//     take:


//     main:
    println "scriptFile: " + workflow.scriptFile
    println "projectDir: " + workflow.projectDir
    println "launchDir: " + workflow.launchDir
    println "workDir: " + workflow.workDir
    println "configFiles: " + workflow.configFiles

        

}


// prepare homoglogous COG group to pp group  mapping in name unsorted manner. so we can know single protein mapping relations
// and here we use homologou pp under eggenog level  max level 2 (root level ) 
// so information stored is indepedent on pps in query species benchmark datasete in which eggno level  could be different 
process homologousPPDetection_COG2PPMapping {
    
    publishDir "${params.homologous_ppPath}", mode: "copy"
    
    
    label "simple_py_process"
    
    
    input: 
        path eggNOG_folder
        path EggNOG_groupPath
        
    output:

    script:
        
    """
        homologous_COG2PP_path="COG2PP/" 

        
        export PYTHONPATH="${projectDir}/../src/utilities/" # "\${projectDir}/../" not working 
        python ${projectDir}/python_scripts/homologousPPDetection_COG2PPMapping.py -s \
        -egg "${eggNOG_folder}/groups/" -t "\${homologous_COG2PP_path}"
    """
    
}

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('-s','--spe_list', type=str, help='spe_list')
    parser.add_argument('-egg','--EggNOG_groupPath', type=str, help='EggNOG_groupPath')
    parser.add_argument('-t','--homologous_COG2PP_path', type=str, help='homologous_COG2PP_path')
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')

eggNOG_group_folder="${eggNOG_folder}/groups"


workflow {
    query2subject_homologousPPDetection_workflow()
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
nextflow run Query_coevolutionComputation_workflow.nf -params-file wc-params.json -c nextflow.config -resume
with "-resume -with-report -with-trace -with-timeline -with-dag dag.png" get more job running report
to view nextflow log file ,  run "ls -lhtra" ,
open file ".nextflow.log"

*/

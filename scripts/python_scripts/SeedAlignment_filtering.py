import argparse
import os 
import glob
import multiprocessing as mp
from Bio import SeqIO
from os import listdir
from os.path import isfile, join
import pandas as pd 
import re

def getSpe2pro2seq_dict(record):
    origSTRINGBacteriaProSeqPath,SpeName=record
    Spe_index=SeqIO.index(os.path.join(origSTRINGBacteriaProSeqPath,SpeName+".fa"),"fasta")
    Spe_index_dict=dict(Spe_index)
    Spe_index.close()
    return((SpeName,Spe_index_dict))


# keep this funcion heree as it need very larget global variable 
# and using subprocess seem not very fast 
def newSingleMSA_filterPhmmer_EggNOG_OrthologousGroup_faidx(file_name):
    try:
        reTemp=re.compile("([0-9]*)\.*")
        init_phmmer_result=pd.read_csv(os.path.join(currentSpe_phmmer_outPath,file_name), 
                                       comment="#", delim_whitespace=True,header=None,index_col=None,engine="python")
        init_phmmer_result['qcov']=(init_phmmer_result.iloc[:,16]-init_phmmer_result.iloc[:,15]+1)/init_phmmer_result.iloc[:,5] # query coverage
        init_OGLen=init_phmmer_result.shape[0]
        # targetname accession   tlen queryname       0-3
        #accession   qlen   E-value  score  bias   4-8
        # #  of  c-Evalue  i-Evalue  score  bias  9-14
        # from    to  from    to    #15-18
        #from    to  acc  description oftarget #19-22
        # checck Tabular output formats in hmmer documentation 
        # here acc is not identitiy .but a similar one
        phmmer_result=init_phmmer_result.loc[init_phmmer_result.iloc[:,21]>0.55,:].copy() # acc > 0.55 acc >0.4
        phmmer_result=phmmer_result.loc[phmmer_result['qcov']>0.8,:] # query coverage > 0.8 >0.65
        if phmmer_result.shape[0]>2500 or phmmer_result.shape[0]>init_OGLen*0.25:
            pass
        else:
            #print("second layer")
            phmmer_result=init_phmmer_result.loc[init_phmmer_result.iloc[:,21]>0.4,:].copy() # acc > 0.55 acc >0.4
            phmmer_result=phmmer_result.loc[phmmer_result['qcov']>0.65,:] # query coverage > 0.8 >0.65
            if phmmer_result.shape[0]>2500 or phmmer_result.shape[0]>init_OGLen*0.25:
                pass
            else:
                #print("third layer")
                phmmer_result=init_phmmer_result.loc[init_phmmer_result.iloc[:,21]>0.25,:].copy() # acc > 0.55 acc >0.4
                phmmer_result=phmmer_result.loc[phmmer_result['qcov']>0.5,:] # query coverage > 0.8 >0.65

        phmmer_result=phmmer_result.loc[phmmer_result.iloc[:,0]!=file_name[0:-10],:] 
    except:
        print("empty file : "+file_name)
        return
        
        #print(file_name)



    #if (phmmer_result.shape[0]>0):
    if (phmmer_result.shape[0]>2500 or phmmer_result.shape[0]>init_OGLen*0.25):
        currentSpe_protein=phmmer_result.iloc[0,3]
        ecoli_dict=SeqIO.index(currentSpeProSeqPath_ByProteins+currentSpe_protein+".fa","fasta")
        ecoli_str=str(ecoli_dict[currentSpe_protein].seq)
        output_str=">"+currentSpe_protein+"\n"+ecoli_str+"\n"

        phmmer_ecoli=phmmer_result.iloc[0,3]
        #print(phmmer_ecoli)
        #???!!! here there are duplicete target proteins as it has  more then two domains aligned 
        # but here we only use unique targert prorteins name and check if duplicted later  
        phmmer_OG_Pros=list(set(phmmer_result.iloc[:,0]))


    # Here is a little tricky
    #which is suppose to slover the problem mention in page 21 of suplemntary files of Qians paperif 
    #First we need to consider alignment(target) overlap in query,  if so, only take on with high identity
    #Second we need to condiser alignent(target) overlap in target, if so. useing intersection of them.
    #or orignal sequence of them. if so. the target sequececn len will increase!!
    #how about only conside the situatin when non overlap in both query and target sequence ?
        for hit_seq in phmmer_OG_Pros: # 
            temp_phmmer_result=phmmer_result.loc[phmmer_result.iloc[:,0]==hit_seq,:].copy()
            if temp_phmmer_result.shape[0]>1:
                temp_phmmer_result.sort_values(by=[15],ascending=True,inplace=True)# make sure segment in query in correct order
                current_index=0
                keepedIdx=list()
                for i in range(1,temp_phmmer_result.shape[0]): # get non overlaped  and order segment in query 
            #         if i>range(temp_phmmer_result.shape[0]):
            #             break
            #        else:
                    if temp_phmmer_result.iloc[current_index,15]<=temp_phmmer_result.iloc[i,15]<=temp_phmmer_result.iloc[current_index,16]:# overlap 
                        #print("Im here")
                        if temp_phmmer_result.iloc[current_index,21]<temp_phmmer_result.iloc[i,21]:
                            #print("Im there")
                            current_index=i
                        if i==(temp_phmmer_result.shape[0]-1) or (temp_phmmer_result.iloc[i+1,15]>temp_phmmer_result.iloc[current_index,16]):
                            keepedIdx.append(current_index)
                    elif temp_phmmer_result.iloc[i,15]>temp_phmmer_result.iloc[current_index,16]:   # none overlap and next segment is after current segemtn 
                  #      keepedIdx.append(current_index) # thi can be omited 
                        keepedIdx.append(i)
                        current_index=i

                keepedIdx=sorted(list(set(keepedIdx)))
                temp_temp_phmmer_result=temp_phmmer_result.iloc[keepedIdx,:].copy()


                current_index2=0
                keepedIdx2=[0]
                for i in range(1,temp_temp_phmmer_result.shape[0]): # get non overlaped and order segments in target
    #                 if i>range(temp_temp_phmmer_result.shape[0]):
    #                     break
    #                 else:
                    if temp_temp_phmmer_result.iloc[current_index,18]<temp_temp_phmmer_result.iloc[i,17]:
                        current_index2=i
                        keepedIdx2.append(current_index2)
                final_temp_phmmer_result=temp_temp_phmmer_result.iloc[keepedIdx2,:]
            else: 
                final_temp_phmmer_result=temp_phmmer_result

            hit_str=str()
            if final_temp_phmmer_result.shape[0]>1:
                print("im non overlaped in both query and target: "+file_name)
                #print(file_name)
            for j in range(final_temp_phmmer_result.shape[0]):
                hit_seq,hit_start,hit_end=final_temp_phmmer_result.iloc[j,[0,17,18]]
                hit_id=reTemp.match(hit_seq).group(1)
                record_dict = species_indexDict_dict[hit_id]
                hit_str=hit_str+str(record_dict[hit_seq].seq)[(hit_start-1):(hit_end-1)]

            output_str=output_str+">"+hit_seq+"\n"+hit_str+"\n"
        out_str_file=open(currentSpe_phmmer_OrthologousGroup_path+currentSpe_protein+".fa","w")
        out_str_file.write(output_str)
        out_str_file.close()

    else:
        print("not enough seq in phmmer seq alignment ",file_name+":",phmmer_result.shape)

if __name__ == '__main__':
    # then filter phmmer resutls using next two  blocks
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-b','--currentSpeProSeqPath_ByProteins', type=str, help='currentSpeProSeqPath_ByProteins')
    parser.add_argument('-pf','--currentSpe_phmmer_OrthologousGroup_path', type=str, help='currentSpe_phmmer_OrthologousGroup_path')
    parser.add_argument('-o','--origSTRINGBacteriaProSeqPath', type=str, help='origSTRINGBacteriaProSeqPath')
    parser.add_argument('-p','--currentSpe_phmmer_outPath', type=str, help='currentSpe_phmmer_outPath') 
    parser.add_argument('-n','--mp_task_nums', type=str, help='mp_task_nums')

    args = parser.parse_args()
    currentSpeProSeqPath_ByProteins=args.currentSpeProSeqPath_ByProteins
    currentSpe_phmmer_OrthologousGroup_path=args.currentSpe_phmmer_OrthologousGroup_path
    origSTRINGBacteriaProSeqPath=args.origSTRINGBacteriaProSeqPath
    currentSpe_phmmer_outPath=args.currentSpe_phmmer_outPath
    mp_task_nums=int(args.mp_task_nums)

    species_fa=glob.glob(os.path.join(origSTRINGBacteriaProSeqPath,"*.fa"))
    species_names=[os.path.basename(f) for f in species_fa]
    species_names=[f[:-3] for f in species_names]
    print("len(species_names):",len(species_names))




    species_names_ArgForGetSpe2pro2seq=[(origSTRINGBacteriaProSeqPath,SpeName) for SpeName in species_names]
    pool = mp.Pool(mp_task_nums)
    getSpe2pro2seq_dict_results=pool.map(getSpe2pro2seq_dict, species_names_ArgForGetSpe2pro2seq)
    pool.close()
    species_indexDict_dict=dict(getSpe2pro2seq_dict_results)
    # CPU times: user 4min 58s, sys: 30.2 s, total: 5min 28s
    # Wall time: 5min 36s

    #run once 
    # now for one species around 3-4 minutes
    # here many operation involde with find one particular protei seq from on species, 
    # better to read them in one run 
    # after read all spefice fasta file, one run takes  3.81 seconds 
    currentSpe_phmmer_files = [f for f in listdir(currentSpe_phmmer_outPath) if isfile(join(currentSpe_phmmer_outPath, f))]
    currentSpe_phmmer_files=[f for f in currentSpe_phmmer_files if "domtblout" in f]

    print("len(currentSpe_phmmer_files):",len(currentSpe_phmmer_files))
    pool = mp.Pool(mp_task_nums) # here dont use many process here as following functin require large memory 
    pool.map(newSingleMSA_filterPhmmer_EggNOG_OrthologousGroup_faidx, currentSpe_phmmer_files)
    pool.close()

    # CPU times: user 5min 45s, sys: 37.8 s, total: 6min 23s
    # Wall time: 13min 46s


# need enviroment py37_pydca

import os
import glob
import sys
import pickle
import time
from datetime import datetime
import pandas as pd

import multiprocessing as mp
from multiprocessing import get_context

from pydca.meanfield_dca import meanfield_dca

if True:
    sys.path.append("../src/utilities")

    from collect_topCoEvos import get_maxBetValue_dict_pydcaFNAPC_array_npz
    from DCA_computation import pydca_mfdca_FN_compresse


os.environ["OMP_NUM_THREADS"] = "4"  # export OMP_NUM_THREADS=4

os.environ["MKL_NUM_THREADS"] = "4"  # export MKL_NUM_THREADS=6

os.environ["NUMBA_NUM_THREADS"] = "4"  # export NUMBA_NUM_THREADS=6


# "/mnt/mnemo6/tao/", /mnt/mnemo6/tao/notebook_data/
notebookData_folder = "/mnt/mnemo6/tao/"
newSTRING_rootFolder = f"{notebookData_folder}PPI_Coevolution/STRING_data_11.5/"
CoEvo_data_folder = f"{notebookData_folder}PPI_Coevolution/CoEvo_data_STRING11.5/"

leftPhylum_tuple = ('1224', '511145')
leftPhylum_EggNOG_maxLevel, leftPhylum_currentSpe_TaxID = leftPhylum_tuple
leftPhylum_input_root_folder = CoEvo_data_folder+leftPhylum_currentSpe_TaxID + \
    "_EggNOGmaxLevel"+leftPhylum_EggNOG_maxLevel+"_eggNOGfilteredData/"
leftPhylum_currentSpeMiddleDataPath = newSTRING_rootFolder+leftPhylum_currentSpe_TaxID + \
    "_EggNOGmaxLevel"+leftPhylum_EggNOG_maxLevel+"_MiddleData/"
print(leftPhylum_currentSpeMiddleDataPath)
with open(leftPhylum_currentSpeMiddleDataPath+'fasta_protein_lens_dict.pickle', 'rb') as handle:
    leftPhylum_fasta_protein_lens = pickle.load(handle)


allPhylum_tuple = ('2', '511145')

allPhylum_EggNOG_maxLevel, allPhylum_currentSpe_TaxID = allPhylum_tuple
allPhylum_input_root_folder = CoEvo_data_folder+allPhylum_currentSpe_TaxID + \
    "_EggNOGmaxLevel"+allPhylum_EggNOG_maxLevel+"_eggNOGfilteredData/"
allPhylum_currentSpeMiddleDataPath = newSTRING_rootFolder+allPhylum_currentSpe_TaxID + \
    "_EggNOGmaxLevel"+allPhylum_EggNOG_maxLevel+"_MiddleData/"
print(allPhylum_currentSpeMiddleDataPath)
with open(allPhylum_currentSpeMiddleDataPath+'fasta_protein_lens_dict.pickle', 'rb') as handle:
    allPhylum_fasta_protein_lens = pickle.load(handle)

allPhylumRandomised_SamePosandNeg_pairedMSA_Nf90_folder = leftPhylum_input_root_folder + \
    "allPhylumRandomised_pair_MSA_Nf90_PasteAlign/"
print(allPhylumRandomised_SamePosandNeg_pairedMSA_Nf90_folder)

DCA_coevolutoin_path = leftPhylum_input_root_folder + \
    "allPhylumRandomised_coevolutoin_result_DCA/"
MI_coevolutoin_path = leftPhylum_input_root_folder + \
    "allPhylumRandomised_coevolutoin_result_MI/"

Benchmark_folder = leftPhylum_input_root_folder+"allPhylumRandomised_Benchmark/"

print(DCA_coevolutoin_path,
      MI_coevolutoin_path, Benchmark_folder)

if not os.path.exists(DCA_coevolutoin_path):
    os.makedirs(DCA_coevolutoin_path)

if not os.path.exists(MI_coevolutoin_path):
    os.makedirs(MI_coevolutoin_path)

if not os.path.exists(Benchmark_folder):
    os.makedirs(Benchmark_folder)


input_msaFiles = glob.glob(
    allPhylumRandomised_SamePosandNeg_pairedMSA_Nf90_folder+"*.fasta")
print(len(input_msaFiles))
input_msaFiles = [os.path.basename(f)[0:-len(".fasta")]
                  for f in input_msaFiles]
input_PPs = [f.split("and") for f in input_msaFiles]


existed_pydcaFNAPC_files = glob.glob(
    DCA_coevolutoin_path+"*_pydcaFNAPC_array.npz")
print("len(existed_pydcaFNAPC_files):", len(existed_pydcaFNAPC_files))


existed_pydcaFNAPC_files = [os.path.basename(
    f) for f in existed_pydcaFNAPC_files]
existed_pydcaFNAPC_pps = [f.split("and") for f in existed_pydcaFNAPC_files]
existed_pydcaFNAPC_pps = [(p1, p2[:-21])
                          for p1, p2 in existed_pydcaFNAPC_pps]
existed_pydcaFNAPC_pp_dict = dict(
    [(pp, 1) for pp in existed_pydcaFNAPC_pps])


ArgForDCA = [
    (p1, p2, leftPhylum_fasta_protein_lens[p1], allPhylum_fasta_protein_lens[p2], allPhylumRandomised_SamePosandNeg_pairedMSA_Nf90_folder, DCA_coevolutoin_path) for p1, p2 in input_PPs if (p1, p2) not in existed_pydcaFNAPC_pp_dict]
print("len(ArgForDCA):", len(ArgForDCA))





if len(ArgForDCA) > 0:

    now = datetime.now()
    print("DCA comp stat now =", now)
    pool = mp.Pool(30)
    dca_results = pool.map(pydca_mfdca_FN_compresse, ArgForDCA)

    time.sleep(10)
    pool.close()

now = datetime.now()
print("DCA comp end now =", now)

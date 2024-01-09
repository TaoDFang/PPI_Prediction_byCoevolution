from Bio import AlignIO
import random
from scipy.spatial.distance import pdist, squareform
from scipy.stats import entropy
from collections import defaultdict
import numpy as np

from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO

import matplotlib.pyplot as plt


################
# note: if you are modifying the alphabet
# make sure last character is "-" (gap)
################
alphabet = "ARNDCQEGHILKMFPSTWYV-"
states = len(alphabet)
a2n = {}
for a, n in zip(alphabet, range(states)):
    a2n[a] = n
################


def aa2num(aa):
    '''convert aa into num'''
    if aa in a2n:
        return a2n[aa]
    else:
        return a2n['-']

# from fasta


def parse_fasta(filename, limit=-1):
    '''function to parse fasta'''
    header = []
    sequence = []
    lines = open(filename, "r")
    for line in lines:
        line = line.rstrip()
        if line[0] == ">":
            if len(header) == limit:
                break
            header.append(line[1:])
            sequence.append([])
        else:
            sequence[-1].append(line)
    lines.close()
    sequence = [''.join(seq) for seq in sequence]
    return np.array(header), np.array(sequence)


def get_dist(msa):
    msa_dists = squareform(pdist(msa, "hamming"))
    return (msa_dists)


def get_msa_colentropy(msa):
    msa_fre = np.empty((states, msa.shape[1]))
    msa_fre[:] = np.nan
    for j in range(msa.shape[1]):
        msa_fre[:, j] = np.bincount(msa[:, j], minlength=states)
    msa_colentropy = entropy(msa_fre, base=states, axis=0)

    return (msa_colentropy)


def get_msa_statistics_dict(l):
    '''converts list of sequences to msa'''
    pp1, pp2, msa_path = l
    msa_filename = msa_path+pp1+"and"+pp2+".fasta"
    names, seqs = parse_fasta(msa_filename)

    msa_ori = []
    for seq in seqs:
        msa_ori.append([aa2num(aa) for aa in seq])
    msa_ori = np.array(msa_ori)

    # compute seq distancre
    msa_dists = get_dist(msa_ori)

    # compute colentropy
    msa_colentropy = get_msa_colentropy(msa_ori)

    mean_msa_dists = np.mean(np.triu(msa_dists, 0))
    mean_msa_colentropy = np.mean(msa_colentropy)

    msa_colGapPercent = np.sum(msa_ori == 20, axis=0)/msa_ori.shape[0]
    mean_msa_colgappercent = np.mean(msa_colGapPercent)

    return (pp1, pp2, msa_ori.shape[0], msa_ori.shape[1], mean_msa_dists, mean_msa_colentropy, mean_msa_colgappercent)


def downsample_msa(l):
    pp1, pp2, inputmsa_path, outputmsa_path, base_msanum = l
    msa = AlignIO.read(inputmsa_path+pp1+"and"+pp2+".fasta", "fasta")
    random.seed(0)

    if len(msa) < base_msanum:
        downsample_idx = random.sample(range(len(msa)), len(msa))
    else:
        downsample_idx = random.sample(range(len(msa)), base_msanum)
    downsample_msa = MultipleSeqAlignment([])

    if base_msanum == 1:
        idx = downsample_idx[0]
        downsample_msa.append(msa[idx, :])
        downsample_msa.append(msa[idx, :])
        downsample_msa.append(msa[idx, :])
    else:
        for idx in downsample_idx:
            downsample_msa.append(msa[idx, :])

    AlignIO.write(downsample_msa, outputmsa_path +
                  pp1+"and"+pp2+".fasta", "fasta")


def onePhylumRandomised_msa(l):
    pp1, pp2, len1, inputmsa_path, outputmsa_path = l

    try:
        msa = AlignIO.read(inputmsa_path+pp1+"and"+pp2+".fasta", "fasta")
    except:
        # this protein pari might not be able to pass filtering step later
        return None

    random.seed(0)
    downsample_idx = random.sample(range(len(msa)), len(msa))

    randomise_msa = MultipleSeqAlignment([])

    for i in range(0, len(msa)):
        idx = downsample_idx[i]
        randomise_msa.append(msa[i, 0:len1]+msa[idx, len1:])

    AlignIO.write(randomise_msa, outputmsa_path +
                  pp1+"and"+pp2+".fasta", "fasta")
    # return (msa, randomise_msa)


def allPhylumRandomised_msa(l):
    pp1, pp2, leftPhylum_len1, allPhylum_len1, leftPhylum_inputmsa_path, allPhylum_inputmsa_path,  outputmsa_path = l

    try:
        leftPhylum_msa = AlignIO.read(
            leftPhylum_inputmsa_path+pp1+"and"+pp2+".fasta", "fasta")
    except:
        # this protein pari might not be able to pass filtering step later
        # print("leftPhylum_error")
        return None

    try:
        allPhylum_msa = AlignIO.read(
            allPhylum_inputmsa_path+pp1+"and"+pp2+".fasta", "fasta")
    except:
        # this protein pari might not be able to pass filtering step later
        print("allPhylum_error")  # check about many is caused by all level
        return None

    random.seed(0)
    if len(leftPhylum_msa) <= len(allPhylum_msa):
        downsample_idx = random.sample(
            range(len(allPhylum_msa)), len(leftPhylum_msa))
    else:
        leftPhylum_msa = leftPhylum_msa[0:len(allPhylum_msa), :]
        downsample_idx = random.sample(
            range(len(allPhylum_msa)), len(leftPhylum_msa))

    randomise_msa = MultipleSeqAlignment([])

    for i in range(0, len(leftPhylum_msa)):
        idx = downsample_idx[i]
        randomise_msa.append(
            leftPhylum_msa[i, 0:leftPhylum_len1]+allPhylum_msa[idx, allPhylum_len1:])

    AlignIO.write(randomise_msa, outputmsa_path +
                  pp1+"and"+pp2+".fasta", "fasta")

    # return (leftPhylum_msa,allPhylum_msa, randomise_msa)

from Bio import SeqIO
import gzip
import pickle
import numpy as np
from collections import defaultdict
import os.path
from datetime import datetime
import functools
import operator

now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


def create_map(data, save=True, k=12):
    path = "map_" + str(k) + ".p"
    if os.path.isfile(path):
        print("exists.")
        with open(path, 'rb') as fp:
            genome_map = pickle.load(fp)
        return genome_map
    # nums = set()
    # for reads in [reads1, reads2]:
    #     for s in reads:
    #         nums.add(len(s))
    genome_map = defaultdict(set)

    for i in range(len(data)):
        if i % 100000 == 0:
            print((i + 1) / len(data), flush=True)
        genome_map[data[i:i + k]].add(i)
    if save:
        with open('map_' + str(k) + '.p', 'wb') as fp:
            pickle.dump(genome_map, fp, protocol=pickle.HIGHEST_PROTOCOL)
    return genome_map


def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)


def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = reverse(seq)
    seq = complement(seq)
    return seq


def get_reads():
    reads1 = []
    reads2 = []
    with gzip.open('ATAC.chr19.R01.fastq.gz', 'rb') as read1:
        for seq1 in read1:
            seq1 = next(read1)
            reads1.append(str(seq1)[2:-3].upper())
            next(read1)
            next(read1)
    with gzip.open('ATAC.chr19.R02.fastq.gz', 'rb') as read2:
        for seq2 in read2:
            seq2 = next(read2)
            temp_seq = str(seq2)[2:-3].upper()
            # reads2.append(temp_seq)
            reads2.append(temp_seq)
            next(read2)
            next(read2)
    return reads1, reads2


def get_data(file):
    for seq_record_chr19 in SeqIO.parse(file, "fasta"):
        print(seq_record_chr19.id, flush=True)
        print(repr(seq_record_chr19.seq), flush=True)
        print(len(seq_record_chr19), flush=True)
    data = seq_record_chr19.seq._data.upper()
    return data


ch19 = "hg19-chr19.fa"
print("getting data...")
data = get_data(ch19)
print("getting reads...")
reads1, reads2 = get_reads()
print("getting map...")

k = 16

genome_map = create_map(data, k=k)
pos = np.full(len(reads1), -1)
for j, read in enumerate(reads1):
    l = []
    # print(j)
    for i in range(len(read) - k):
        l.append(list(np.array(list(genome_map[read[i:i + k]])) - i))
    l2 = functools.reduce(operator.iconcat, l, [])
    if len(l2) == 0:
        continue
    # print()
    # print(l2)
    idx = np.argmax(np.bincount(np.array(l2)))
    pos[j] = idx
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
del genome_map

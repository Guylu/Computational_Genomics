from Bio import SeqIO
import gzip
import pickle
import numpy as np
from collections import defaultdict
import os.path
from datetime import datetime
import time
import argparse
from matplotlib import pyplot as plt


def create_map(data, save=True, k=12):
    path = "map_" + str(k) + ".p"
    if os.path.isfile(path):
        print("exists.")
        with open(path, 'rb') as fp:
            genome_map = pickle.load(fp)
        return genome_map
    genome_map = defaultdict(set)

    for i in range(len(data)):
        genome_map[data[i:i + k]].add(i)
    if save:
        with open('map_' + str(k) + '.p', 'wb') as fp:
            pickle.dump(genome_map, fp, protocol=pickle.HIGHEST_PROTOCOL)
    return genome_map


complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}


def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = list(seq[::-1])
    return ''.join([complement_dict[base] for base in seq])


def get_reads(path, rv=False):
    reads1 = []
    if rv:
        with gzip.open(path, 'rb') as read1:
            for seq1 in read1:
                seq1 = next(read1)
                t = str(seq1)[2:-3].upper()
                reads1.append(reverse_complement(t))
                next(read1)
                next(read1)
    else:
        with gzip.open(path, 'rb') as read1:
            for seq1 in read1:
                seq1 = next(read1)
                reads1.append(str(seq1)[2:-3].upper())
                next(read1)
                next(read1)
    return np.array(reads1)


def get_data(file):
    for seq_record_chr19 in SeqIO.parse(file, "fasta"):
        pass
    data = seq_record_chr19.seq._data.upper()
    for i in range(len(data)):
        if data[i] != 'N':
            break
    for j in range(1, len(data)):
        if data[-j] != 'N':
            break
    print(i, -j + 1)
    global cut_start, cut_end
    cut_start = i
    cut_end = -j + 1
    return data[i:-j + 1]


def vec_translate(a, my_dict):
    return np.vectorize(my_dict.__getitem__, otypes=[np.object])(a)


def vec_translate_offset(a, my_dict, k):
    d = lambda key: set(np.array(list(my_dict[key])) - k)
    return np.vectorize(d, otypes=[np.object])(a)


def slicer_vectorized(a, start, end):
    b = a.view((str, 1)).reshape(len(a), -1)[:, start:end]
    return np.fromstring(b.tostring(), dtype=(str, end - start))


def map_me_up(genome_map, data, k, path):
    reads1 = get_reads(path)
    reads1rc = get_reads(path, True)
    # print("mapping...")
    start2 = time.time()
    read1_len = len(reads1[0])

    maps11 = vec_translate(slicer_vectorized(reads1, 0, k), genome_map)
    maps12 = vec_translate_offset(slicer_vectorized(reads1, read1_len - k, read1_len), genome_map, read1_len - k)

    union_sets = np.vectorize(set.union)
    unionised = union_sets(maps11, maps12)

    del maps11, maps12

    #######

    maps11rc = vec_translate(slicer_vectorized(reads1rc, 0, k), genome_map)
    maps12rc = vec_translate_offset(slicer_vectorized(reads1rc, read1_len - k, read1_len), genome_map,
                                    read1_len - k)

    unionisedrc = union_sets(maps11rc, maps12rc)
    del maps11rc, maps12rc

    for j, idx in enumerate(unionised):
        if len(idx) == 2:
            idx = list(idx)
            score1 = sum(a == b for a, b in zip(data[idx[0]:idx[0] + read1_len], reads1[j]))
            score2 = sum(a == b for a, b in zip(data[idx[1]:idx[1] + read1_len], reads1[j]))
            if score1 < score2:
                unionised[j].remove(idx[0])
            else:
                unionised[j].remove(idx[1])

    for j, idx in enumerate(unionisedrc):
        if len(idx) == 2:
            idx = list(idx)
            score1 = sum(a == b for a, b in zip(data[idx[0]:idx[0] + read1_len], reads1rc[j]))
            score2 = sum(a == b for a, b in zip(data[idx[1]:idx[1] + read1_len], reads1rc[j]))
            if score1 < score2:
                unionisedrc[j].remove(idx[0])
            else:
                unionisedrc[j].remove(idx[1])

    for index, (u, urc) in enumerate(zip(unionised, unionisedrc)):
        if len(u) == 0 and len(urc) == 1:
            unionised[index] = urc
        elif len(u) == len(urc) == 1:
            i1, i2 = list(u)[0], list(urc)[0]
            score1 = sum(a == b for a, b in zip(data[i1:i1 + read1_len], reads1[index]))
            score2 = sum(a == b for a, b in zip(data[i2:i2 + read1_len], reads1rc[index]))
            if score1 < score2:
                unionised[index].remove(i1)
                unionised[index].add(i2)
    return unionised


parser = argparse.ArgumentParser()
parser.add_argument('k', help='k-mer', default=20)
command_args = parser.parse_args()
k = int(command_args.k)
print(k)
ch19 = "hg19-chr19.fa"
# print("getting data...")
data = get_data(ch19)  # [100000:200000]
# print("getting map...")
t1 = time.time()
genome_map = create_map(data, k=k)
print("Generating Dict Time:" + str(time.time() - t1))
# print("getting reads...")
gene_maps = []
for path in ['ATAC.chr19.R01.fastq.gz', 'ATAC.chr19.R02.fastq.gz']:
    print(path + ":")
    start = time.time()
    gene_maps.append(map_me_up(genome_map, data, k, path))
    end = time.time()
    print("Algo Time =", end - start)
    print()
    m = 0
    e = 0
    d = 0
    for l in gene_maps[-1]:
        if len(l) == 1:
            m += 1
        elif len(l) == 0:
            e += 1
        else:
            d += 1
    # print(m, e, d)
    print("Success Rate: ")
    print(m / (m + e + d))
    print()
    print()
    print()
np.save('uni' + str(k) + '_1.npy', gene_maps[0])
np.save('uni' + str(k) + '_2.npy', gene_maps[1])
hist = np.zeros(len(data) + cut_start - cut_end)
for m in gene_maps:
    for s in m:
        if len(s) == 1:
            hist[list(s)[0] - cut_start] += 1

plt.figure(figsize=(40, 10))
plt.plot(hist)
plt.title("histogram of reads on CHR19 with k=" + str(k))
plt.savefig('graphs of ' + str(k) + '.png')
plt.savefig('graphs of ' + str(k) + '_test.png', bbox_inches='tight')
plt.show()
del genome_map

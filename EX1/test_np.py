from Bio import SeqIO
import gzip
import pickle
import numpy as np
from collections import defaultdict
import os.path


def create_map(data, save=True, k=12):
    path = "map_" + str(k) + ".p"
    if os.path.isfile(path):
        print("exists.")
        with open(path, 'rb') as fp:
            genome_map = pickle.load(fp)
        print("returning map")
        return genome_map
    # nums = set()
    # for reads in [reads1, reads2]:
    #     for s in reads:
    #         nums.add(len(s))
    genome_map = defaultdict(list)

    for i in range(len(data)):
        if i % 50000 == 0:
            print((i + 1) / len(data), flush=True)
        genome_map[data[i:i + k]].append(i)
    for k in genome_map.keys():
        genome_map[k] = np.array(genome_map[k])
    if save:
        with open('map_' + str(k) + '.p', 'wb') as fp:
            pickle.dump(genome_map, fp, protocol=pickle.HIGHEST_PROTOCOL)
    print("returning map")
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
            reads2.append(reverse_complement(temp_seq))
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

k = 35

genome_map = create_map(data, k=k)


def vec_translate(a, my_dict):
    return np.vectorize(my_dict.__getitem__, otypes=[np.object])(a)


def vec_ret(read):
    THR = 20
    indices = list(genome_map[read])
    scores = np.zeros(len(indices))
    s = set()
    if len(indices) == 0: return s
    for t, index in enumerate(indices):
        scores[t] = np.count_nonzero(
            np.array(list(read)) == np.array(list(data[index:index + len(read)])))
    sc = np.argmax(scores)
    best_fit_score = scores[sc]
    if best_fit_score > THR:
        best_fit = indices[sc]
        s.add(best_fit)
    return s


def slicer_vectorized(a, start, end):
    b = a.view((str, 1)).reshape(len(a), -1)[:, start:end]
    return np.fromstring(b.tostring(), dtype=(str, end - start))


r1 = np.array(list(reads1))
maps1 = vec_translate(slicer_vectorized(r1, 0, k), genome_map)
m, e, d = 0, 0, 0
for l in maps1:
    if len(l) == 1:
        m += 1
    elif len(l) == 0:
        e += 1
    else:
        d += 1
print(m, e, d)
#
# res = np.full(len(reads1), -1)
# for j, read in enumerate(reads1):
#     if j % 100 == 0:
#         print(j)
#     practical_len = len(read) - k
#     l = []
#     for i in range(practical_len):
#         idxs = genome_map[read[i:i + k]]
#         if len(idxs) > 0:
#             l.append(idxs - i)
#     if len(l) / practical_len > 0.7:
#         res[j] = np.argmax(np.bincount(np.concatenate(l)))
#
del genome_map

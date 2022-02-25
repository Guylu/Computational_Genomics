from Bio import SeqIO
import gzip
import pickle
import numpy as np
from collections import defaultdict
import os.path
from datetime import datetime

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
        if i % 50000 == 0:
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

k = 20

genome_map = create_map(data, k=k)


# reads1, reads2 = reads1[:100], reads2[:100]


# genome_map = create_map(data, k=5)

##### K MERS:

#
# print("mapping...")
# pos = [""] * len(data)
# yes1, yes2, no1, no2 = 0, 0, 0, 0
# THR = 20
# for i, reads in enumerate([reads1, reads2]):
#     for j, read in enumerate(reads):
#         # print("reads:" + str(i))
#         # print(j)
#         if j % 100 == 0:
#             print(j)
#         for l in range(len(read) - k):  # +1?
#             try:
#                 indices = genome_map[read[l:l + k]]
#                 scores = np.zeros(len(indices))
#                 # r = np.array(list(read))
#                 # t = np.tile(r,(len(indices),1))
#                 # data = np.array(list(data))
#                 # indices = np.array(indices)
#                 # scores = np.count_nonzero(t == np.array(data[indices:indices + len(read) + 1]))
#                 for t, index in enumerate(indices):
#                     scores[t] = np.count_nonzero(
#                         np.array(list(read)) == np.array(list(data[index:index + len(read)])))
#                 best_fit_score = np.max(scores)
#                 if best_fit_score < THR:
#                     continue
#                 best_fit = indices[np.argmax(scores)]
#                 pos[best_fit + l] = read
#                 if i == 0:
#                     yes1 += 1
#                 else:
#                     yes2 += 1
#             except KeyError as e:
#                 # print(read + " not found")
#                 # if i == 0:
#                 #     no1 += 1
#                 # else:
#                 #     no2 += 1
#                 pass
# print(yes1, flush=True)
# print(yes2, flush=True)
#

##### 20 mapping

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


print("mapping...")
pos = [""] * len(data)
yes1, yes2, no1, no2 = 0, 0, 0, 0

r1 = np.array(list(reads1))

read1_len = len(r1[0])


maps11 = vec_translate(slicer_vectorized(r1, 0, k), genome_map)
maps12 = vec_translate(slicer_vectorized(r1, read1_len - k, read1_len), genome_map)

maps12 = np.delete(maps12, [i for i in range(read1_len - k)])
maps12 = np.append(maps12, [set() for _ in range(read1_len - k)])

union_sets = np.vectorize(set.union)
unionised = union_sets(maps11, maps12)

dups = []
for j, idx in enumerate(unionised):
    if len(idx) == 2:
        dups.append([j, list(idx)])
# j is the js reads and idx is the index
for j, idx in dups:
    score1 = sum(a == b for a, b in zip(data[idx[0]:idx[0] + read1_len], reads1[j]))
    score2 = sum(a == b for a, b in zip(data[idx[1]:idx[1] + read1_len], reads1[j]))
    if score1 < score2:
        unionised[j].remove(idx[0])
    else:
        unionised[j].remove(idx[1])

m = 0
e = 0
d = 0
for l in unionised:
    if len(l) == 1:
        m += 1
    elif len(l) == 0:
        e += 1
    else:
        d += 1
print(m, e, d)

# v = np.vectorize(vec_ret)
# test = v(slicer_vectorized(reads1, 0, k))


# for i, reads in enumerate([reads1, reads2]):
#     for j, read in enumerate(reads):
#         index = genome_map[read[: k]]
#         if len(index) == 1:
#             pos[index[0]] = read
#             if i == 0:
#                 yes1 += 1
#             else:
#                 yes2 += 1


del genome_map
# maps2 = vec_translate(reads2, genome_map)

#
# print("mapping...")
# pos = [""] * len(data)
# yes1, yes2, no1, no2 = 0, 0, 0, 0
# for i, reads in enumerate([reads1, reads2]):
#     for j, read in enumerate(reads):
#         # print("reads:" + str(i))
#         # print(j)
#         try:
#             index = genome_map[read[: k]]
#             if len(index) == 1:
#                 pos[index[0]] = read
#                 if i == 0:
#                     yes1 += 1
#                 else:
#                     yes2 += 1
#         except KeyError as e:
#             # print(read + " not found")
#             if i == 0:
#                 no1 += 1
#             else:
#                 no2 += 1
#             pass

# k = 200
# for reads in [reads1, reads2]:
#     for i, seq in enumerate(reads):
#         print((i + 1) / len(reads1))
#         if seq[:k] not in genome_map.keys():
#             genome_map[seq[:k]] = []
#         genome_map[seq[:k]].append(i)


#
#    Basic Algo:
# for i in data:
#     get seq on lenght k at position i
#     if dict[seq] = None:
#       dict[seq] = []
#     dict[seq].append(i)


#   how to retrieve?
# for read in reads:
#   pos_of_read_[i] = dict[read] with the least amount of error

now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#!/usr/bin/env python

import collections, sys
import time
from Bio import Seq, SeqIO, SeqRecord

k = 19

def contig_to_string(c):
    return c[0] + ''.join(x[-1] for x in c[1:])

def get_contig(d,km):
    c_fw = get_contig_forward(d,km)

    c_bw = get_contig_forward(d,twin(km))

    if km in fw(c_fw[-1]):
        c = c_fw
    else:
        c = [twin(x) for x in c_bw[-1:0:-1]] + c_fw
    return contig_to_string(c),c


def get_contig_forward(d,km):
    c_fw = [km]

    while True:
        if sum(x in d for x in fw(c_fw[-1])) != 1:
            break

        cand = [x for x in fw(c_fw[-1]) if x in d][0]
        if cand == km or cand == twin(km):
            break # break out of cycles or mobius contigs
        if cand == twin(c_fw[-1]):
            break # break out of hairpins

        if sum(x in d for x in bw(cand)) != 1:
            break

        c_fw.append(cand)

    return c_fw

def fw(km):
    for x in 'ACGT':
        yield km[1:]+x

def bw(km):
    for x in 'ACGT':
        yield x + km[:-1]

def twin(km):
    return Seq.reverse_complement(km)

def print_dbg(cs):
    """ Print out in Fasta format """
    for i,x in enumerate(cs):
        print('>contig%d\n%s\n'%(i,x))

# def all_contigs(d,k):
done = set()
r = []
d = collections.defaultdict(int)

for line in sys.stdin:
    s_line = line.split('\t')
    d[s_line[0].strip()] = s_line[1].strip()
for x in d:
    if x not in done:
        s,c = get_contig(d,x)
        for y in c:
            done.add(y)
            done.add(twin(y))
        r.append(s)

G = {}
heads = {}
tails = {}
for i,x in enumerate(r):
    G[i] = ([],[])
    heads[x[:k]] = (i,'+')
    tails[twin(x[-k:])] = (i,'-')

for i in G:
    x = r[i]
    for y in fw(x[-k:]):
        if y in heads:
            G[i][0].append(heads[y])
        if y in tails:
            G[i][0].append(tails[y])
    for z in fw(twin(x[:k])):
        if z in heads:
            G[i][1].append(heads[z])
        if z in tails:
            G[i][1].append(tails[z])
print_dbg(r)
    # return G,r

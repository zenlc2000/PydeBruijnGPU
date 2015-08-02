#!/usr/bin/env python

import collections,sys
import os
from Bio import Seq, SeqIO, SeqRecord


k = 19
limit  = 0

d = collections.defaultdict(int)

def kmers(seq,k):       # possibly good candidate
    for i in xrange(len(seq)-k+1):
        yield seq[i:i+k]

def twin(km):
    return Seq.reverse_complement(km)

reads = SeqIO.parse(sys.stdin,'fastq')
for read in reads:
    seq_s = str(read.seq)
    seq_l = seq_s.split('N')
    for seq in seq_l:
        for km in kmers(seq,k):
            d[km] +=1
        seq = twin(seq)
        for km in kmers(seq,k):
            d[km] += 1
for key, val in d.iteritems():
    print str(key),'\t',str(val)

# d1 = [x for x in d if d[x] <= limit]
# for x in d1:
    # del d[x]

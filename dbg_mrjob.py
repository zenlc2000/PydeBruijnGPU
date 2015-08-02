import collections, sys, re
import time
from Bio import Seq, SeqIO, SeqRecord
from mrjob.job import MRJob
from mrjob.step import MRStep
import mrjob


class Assembly(MRJob):

    reads = []
    d = collections.defaultdict(int)
    k = 19

    INPUT_PROTOCOL = mrjob.protocol.RawValueProtocol

    def configure_options(self):
        super(Assembly, self).configure_options()
        self.add_passthrough_option(
            '--k', help="Specify k-value for the job")


    def twin(self,km):
        return Seq.reverse_complement(km)

    def kmers(self,seq,k):
        for i in xrange(len(seq)-k+1):
            yield seq[i:i+k]

    def fw(self,km):
        for x in 'ACGT':
            yield km[1:]+x

    def bw(self,km):
        for x in 'ACGT':
            yield x + km[:-1]


    # def build(fn,k=31,limit=1):
    def build(self, _, data):
        junk = re.match('.*[^ACGT].*', data.strip())
        if not junk:
            # print junk
            self.reads.append(data.strip())
            # print data
        #  self.reads

        # k = 19
        # # for f in fn:
        # # reads = SeqIO.parse(data,'fastq')
        for read in self.reads:
        #     seq_s = str(read.seq)
        #     seq_l = seq_s.split('N')
            # for seq in read:
            kmer = self.kmers(read,self.k)
            for km in kmer:
                self.d[km] +=1
            read = self.twin(read)
            kmer = self.kmers(read,self.k)
            for km in kmer:
                self.d[km] += 1
                        #
        # d1 = [x for x in d if d[x] <= limit]
        # for x in d1:
        #     del d[x]

        # return d
        # print dict(self.d)
        for key, val in self.d.iteritems():
            yield key, val


    def contig_to_string(self,c):
        return c[0] + ''.join(x[-1] for x in c[1:])

    def get_contig(self,d,km):
        c_fw = get_contig_forward(d,km)

        c_bw = get_contig_forward(d,twin(km))

        if km in fw(c_fw[-1]):
            c = c_fw
        else:
            c = [twin(x) for x in c_bw[-1:0:-1]] + c_fw
        return contig_to_string(c),c


    def get_contig_forward(self,d,km):
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

    def all_contigs(self,d,k):
        k = 19
        done = set()
        r = []
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

        return G,r



    def print_GFA(self,G,cs,k):
        """ Print in GFA format """
        print "H  VN:Z:1.0"
        for i,x in enumerate(cs):
            print "S\t%d\t%s\t*"%(i,x)

        for i in G:
            for j,o in G[i][0]:
                print "L\t%d\t+\t%d\t%s\t%dM"%(i,j,o,k-1)
            for j,o in G[i][1]:
                print "L\t%d\t-\t%d\t%s\t%dM"%(i,j,o,k-1)

    def print_dbg(self,cs):
        """ Print out in Fasta format """
        for i,x in enumerate(cs):
            print('>contig%d\n%s\n'%(i,x))

    def init(self):
        # k = int(sys.argv[1])
        k = self.options.k
        self.reads = []
        # JUNK = re.compile('.*[^ACGT].*', flags=0)
        # t0 = time.clock()
        # tw = time.time()
        # d = build(sys.argv[2:],k,1)
        # print time.clock() - t0, "seconds process time"
        # print time.time() - tw, "seconds wall time"
        # G,cs = all_contigs(d,k)
        # # print_GFA(G,cs,k)
        # print_dbg(cs)

    def steps(self):
        return [
            MRStep(mapper=self.build,
                    reducer_init=self.init,
                   reducer=self.all_contigs) #,
        ]

if __name__ == '__main__':
    Assembly.run()

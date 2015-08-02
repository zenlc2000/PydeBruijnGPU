from mrjob.job import MRJob
from mrjob.step import MRStep
import mrjob, re, collections

class Assemble(MRJob):

    reads = []
    d = collections.defaultdict(int)

    def build(self, _, data):
        junk = re.match('.*[^ACGT].*', data.strip())
        if not junk:
            print data
            self.reads.append(data.strip())

        k = 19
        # for f in fn:
        # reads = SeqIO.parse(data,'fastq')
        # for read in self.reads:
        #     seq_s = str(read.seq)
        #     seq_l = seq_s.split('N')
        #     for seq in seq_l:
        #         for km in kmers(seq,k):
        #             d[km] +=1
        #         seq = twin(seq)
        #         for km in kmers(seq,k):
        #             d[km] += 1
                        #
        # d1 = [x for x in d if d[x] <= limit]
        # for x in d1:
        #     del d[x]

        # return d
        # yield d



    def steps(self):
        return [
            MRStep(mapper=self.build)
            ]


if __name__ == '__main__':
    Assemble(args=None).run()

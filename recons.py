from local_para import ReconsPara
import cPickle as pkl

def writeBinary(f, para, content):
    from struct import pack
    f.write(pack(para, content))

def writeInt(f, content):
    writeBinary(f, 'i', content)

def writeFloat(f, content):
    writeBinary(f, 'f', content)

class ReconsturctionData:

    def __init__(self):
        # head
        self.n_particle = None;
        self.n_strand = None;
        self.factor = None;

        # head2
        self.n_group = None;
        self.n_frame = None;
        self.mcx = None;

        # index computed later

        # guide section

    def computeIndices(self):
        self.idx_guide = 12+128+12+12;
        self.idx_frame = self.idx_guide + 4*self.n_group + \
            24 * self.n_group * self.factor + \
            self.n_frame * (4+64+48*self.n_group*self.factor)


if __name__ == "__main__":

    import sys
    exportName = sys.argv[1]
    assert(exportName != "")

    paras = ReconsPara()
    with (open(exportName, 'wb'), open(paras.reference, 'rb')) as (out, anim):

        infos = pkl.load(open(paras.info))
        weights = pkl.load(open(paras.weights, 'rb'))

        # head
        writeInt(out, infos[1])  # particle
        writeInt(out, infos[0])  # strand
        writeInt(out, infos[2])  # factor

        # head2
        writeInt(out, len(weights[0]))
        writeInt(out, )

        anim.close()

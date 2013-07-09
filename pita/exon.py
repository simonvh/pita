class Exon:
    def __init__(self, chrom, start, end, strand):
        if end < start:
             raise ValueError, "exon end < start (%s:%s%s%s)" % (chrom, start, strand, end)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.evidence = []
        self.validated = False
        self.stats = {}

    def __nonzero__(self):
        return 1
    
    def __len__(self):
        return self.end - self.start

    def __repr__(self):
        return "%s:%s-%s" % (self.chrom, self.start, self.end)

    def add_evidence(self, ev):
        if not ev in self.evidence:
            self.evidence.append(ev)

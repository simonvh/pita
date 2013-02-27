class Exon:
	def __init__(self, chr, start, end, strand):
		if end < start:
			 raise ValueError, "exon end < start (%s:%s%s%s)" % (chr, start, strand, end)
		self.chr = chr
		self.start = start
		self.end = end
		self.strand = strand
		self.evidence = []
		self.links = {}
		self.linked = False
		self.validated = False
		self.stats = {}

	def __nonzero__(self):
		return 1
	
	def __len__(self):
		return self.end - self.start

	def __repr__(self):
		return "%s:%s-%s" % (self.chr, self.start, self.end)

	def add_evidence(self, ev):
		if not ev in self.evidence:
			self.evidence.append(ev)

	def add_link(self, exon, ev):
		if exon.chr != self.chr:
			print "%s:%s-%s %s - %s:%s-%s %s" % (self.chr, self.start, self.end, self.strand, exon.chr, exon.start, exon.end, exon.strand)
			raise ValueError, "Different chromosomes!"
		if exon.start <= self.end:
			print "%s:%s-%s %s - %s:%s-%s %s" % (self.chr, self.start, self.end, self.strand, exon.chr, exon.start, exon.end, exon.strand)
			raise ValueError, "exons overlap, or in wrong order"
			
		if exon.strand != self.strand:
			print "%s:%s-%s %s - %s:%s-%s %s" % (self.chr, self.start, self.end, self.strand, exon.chr, exon.start, exon.end, exon.strand)
			raise ValueError, "strands don't match"
		
		self.links.setdefault(exon, [])
		if not ev in self.links[exon]:
			self.links[exon].append(ev)
		exon.linked = True

	def get_linked_exons(self):
		return self.links.keys()

	def linked_chains(self, ev=None):
		bla = []
		for exon,transcripts in self.links.items():
			if not ev or ev in transcripts:
				chains = exon.linked_chains(ev)
				for chain in chains:
					bla.append([self] + chain)
		if not bla:
			bla = [[self]]
		return bla


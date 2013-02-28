from BCBio import GFF
from pita.collection import *
import pprint
import sys

def _gff_type_iterator(feature, ftype):
	if feature.type == ftype:
		yield feature
	else:
		for feature in feature.sub_features:
			for f in _gff_type_iterator(feature, ftype):
				yield f

def read_gff(fname, format="gtf", collection=None, prefix=None):
	#limits = dict(gff_type = ["mRNA", "exon"])
	smap = {"1":"+",1:"+","-1":"-",-1:"-"}
	link = {}
	for rec in GFF.parse(open(fname)):
		chrom = rec.id
		for feature in rec.features:
			for gene in _gff_type_iterator(feature, 'mRNA'):
				for exon in [f for f in gene.sub_features if f.type == 'exon']:
					#print gene.strand, exon.location, gene.id
					link[gene.id] = link.setdefault(gene.id, 0) + 1
					start = int(exon.location.start.position) - 1	
					end = int(exon.location.end.position)
					strand = smap[exon.strand]
					collection.add(chrom, start, end, exon.strand, "bla", transcript=gene.id)	

	for transcript_id, count in link.items():
		if count > 0:
			collection.link_transcript_exons(transcript_id)
	return collection

def read_bed(bedfile, collection=None, prefix=None):
	names = {}
	link = {}
	for line in open(bedfile):
		if not line.startswith("track"):
			vals = line.strip().split("\t")
			# More than one exon
			chromStart = int(vals[1])
			if int(vals[9]) > 1:
				sizes = [int(x) for x in vals[10].split(",")[:-1]]
				starts = [int(x) for x in vals[11].split(",")[:-1]] 
				i = 1
				name = "%s_%s_%s" % (vals[0], vals[3], i)
				while names.has_key(name):
					i += 1
					name = "%s_%s_%s" % (vals[0], vals[3], i)
				names[name] = 1
					
				for start, size in zip(starts, sizes):
					#print 	
					#print vals[0], chromStart + start, chromStart + start + size, vals[5], name
					if size < 0:
						print line
					collection.add(vals[0], chromStart + start, chromStart + start + size, vals[5], "bla", transcript=name)	

	#print link
	for name in names.keys():
		collection.link_transcript_exons(name)
	return collection


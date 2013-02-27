from pita.collection import *

def read_gff(gfffile, format="gtf", collection=None, prefix=None):
	link = {}
	if format == "gtf":
		sep = " "
		names = ["transcript_id", "name", "featureId"]
	elif format == "gff3" or format == "gff":
		sep = "="
		names = ["Name", "pacid"]
	else:
		raise ValueError, "Unknown format %s" % format

	for ln_i,line in enumerate(open(gfffile)):
		if ln_i % 10000 == 0:
			print "%s lines..." % ln_i
		if not line.startswith("#"):
			vals = line.strip().split("\t")
			if vals[2] == "exon":
				params = dict([[y.strip('"') for y in x.strip().split(sep)] for x in vals[8].strip(";").split(";")])

				name = None
				for x in names:
					if params.has_key(x):
						name = params[x]
						if prefix:
							name = prefix + name
			
				if name:
					link[name] = link.setdefault(name, 0) + 1
					start = int(vals[3])
					if format == "gff" or format == "gff3":
						start -= 1
					collection.add(vals[0], start, vals[4], vals[6], "bla", transcript=name)	

	#print link
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


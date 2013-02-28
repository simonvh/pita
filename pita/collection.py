from pita.exon import *
from pita.io import *
import sys

class Collection:
	def __init__(self):
		self.exons = {}
		self.transcript_exons = {}
		self.exon_index = {}

	def add_annotation(self, fname, fformat, prefix=""):
		sys.stderr.write("Reading %s\n" % fname)
		c = Collection()
		if fformat in ["gtf", "gff3"]:
			c = read_gff(fname, fformat, collection=c, prefix=prefix)
		elif fformat == "bed":
			c = read_bed(fname, collection=c, prefix=prefix)
		sys.stderr.write("Done\n")

		sys.stderr.write("Updating models\n")
		ret = self.update_with(c)
		sys.stderr.write("Done\n")

	def delete_left_exons(self, transcript, n):
		#print "Delete %s exons from the left of %s" % (n, transcript)
		exons = self.get_transcript_exons(transcript)
		for i in range(n):
			exon = exons[i]
			if transcript in exon.evidence:
				exon.evidence.remove(transcript)
			else:
				sys.stderr.write("Exon %s has no evidence %s\n" % (exon, transcript))
			self.transcript_exons[transcript] = self.transcript_exons[transcript][1:]
			if self.exons[exon.chr][exon].has_key(transcript):
				del self.exons[exon.chr][exon][transcript]
			else:
				sys.stderr.write("%s no key %s\n" % (exon, transcript))

			if len(exon.evidence) == 0:
				self.delete_exon(exon)
	
	def delete_right_exons(self, transcript, n):
		#print "Delete %s exons from the right of %s" % (n, transcript)
		exons = self.get_transcript_exons(transcript)
		#print len(exons), len(exons) - 1, len(exons) - n - 1
		#print transcript
		#print exons
		#print len(exons) - n - 1, len(exons) - n, n
		if  exons[len(exons) - n - 1].links.has_key(exons[len(exons) - n]):
			try:
				exons[len(exons) - n - 1].links[exons[len(exons) - n]].remove(transcript)
			except:
				sys.stderr.write("ERROR: right exon\n")
				for e in exons:
					sys.stderr.write("%s\n%s\n" % (str(e), str(e.links)))
		else:
			sys.stderr.write("No link for %s\n" % exons[len(exons) - n])
		for i in range(len(exons) - 1, len(exons) - n - 1, -1):
			exon = exons[i]
			sys.stderr.write("Del %s from %s\n" % (exon, transcript))
			if transcript in exon.evidence:
				exon.evidence.remove(transcript)
			else:
				sys.stderr.write("Exon %s has no evidence %s\n" % (exon, transcript))
			#print self.transcript_exons[transcript]
			self.transcript_exons[transcript] = self.transcript_exons[transcript][:-1]
			#print self.transcript_exons[transcript]

			if self.exons[exon.chr][exon].has_key(transcript):
				del self.exons[exon.chr][exon][transcript]
			else:
				sys.stderr.write("%s no key %s\n" % (exon, transcript))
			if len(exon.evidence) == 0:
				self.delete_exon(exon)

	def delete_exon(self, exon):
		#print "Deleting %s" % exon
		for t in self.get_transcripts(exon):
			sys.stderr.write("t: %s\n" % str(t))
			self.transcript_exons[t].remove(exon)
		
		del self.exons[exon.chr][exon]
		del self.exon_index[exon.chr][(exon.start, exon.end, exon.strand)]
		
		exon = None
		

	def add_exon(self, exon, ev=None, transcript="", exon_id=None):
		#print "Adding %s" % exon
		ev = transcript
		self.exon_index.setdefault(exon.chr, {})	
		if not self.exon_index[exon.chr].has_key((exon.start, exon.end, exon.strand)):
			self.exon_index[exon.chr][(exon.start, exon.end, exon.strand)] = exon
		self.exons.setdefault(exon.chr, {}).setdefault(exon, {})[transcript] = 1
	
		if exon_id:
			exon.exon_id = exon_id
		
		if transcript:
			if exon_id:
				for other_exon in self.get_transcript_exons(transcript):
					if abs(exon_id - other_exon.exon_id) == 1:
						if exon.start < other_exon.start:
							self.link_exons(exon, other_exon, transcript)
						else:
							self.link_exons(other_exon, exon, transcript)
			if not exon in self.get_transcript_exons(transcript):
				self.transcript_exons.setdefault(transcript, []).append(exon)
				self.transcript_exons[transcript] = sorted(self.transcript_exons[transcript], key=lambda x:x.start)	
		
		return exon

	def add(self, chr, start, end, strand, ev="", transcript="", exon_id=""):
		e = self.get_exon(chr, int(start), int(end), strand)
		if not e:
			e = Exon(chr, int(start), int(end), strand)
		e.add_evidence(transcript)
		return self.add_exon(e, ev, transcript, exon_id)
		
	def link_exons(self, e1, e2, ev=None):
		try:
			e1.add_link(e2, ev)
		except:
			sys.stderr.write("refusing to link %s and %s, overlap or wrong order" % (e1, e2))

	def link_transcript_exons(self, transcript):
		#print "HOEIEIEO"
		exons = self.get_transcript_exons(transcript)
		exons = sorted(exons, key=lambda x: x.start)
		#print exons
		for i in range(len(exons) - 1):
			try:
				self.link_exons(exons[i], exons[i + 1], transcript)
			except:
				sys.stderr.write("COULD NOT LINK %s and %s\n" % (exons[i], exons[i + 1]))

	def get_exons(self, chr=None):
		if chr:
			if self.exons.has_key(chr):
				return self.exons[chr].keys()
			else:
				return []
		else:
			exons = []
			for es in self.exons.values():
				for e in es.keys():
					#print "e",e
					exons.append(e)
			return exons
	
	
	def get_exon(self, chr, start, end, strand=None):
		if not self.exons.has_key(chr):
			return None

		if self.exon_index[chr].has_key((start, end, strand)):
			return 	self.exon_index[chr][(start, end, strand)]
		
		#for coords, exon in self.exon_index[chr].items():
		#	#if strand:
		#	#	if coords == (start, end, strand):
		#	#		return exon
		#	#else:
		#	if 1:
		#		if coords[0] == start and coords[1] == end:
		#			if not strand:
		#				return exon
		#			elif coords[2] == strand:
		#				return exon
	
	def get_transcript_exons(self, transcript):
		if self.transcript_exons.has_key(transcript):
			return sorted(self.transcript_exons[transcript], key=lambda x: x.start)
		return []

	def get_transcript_start_exon(self, transcript):
		exons = self.get_transcript_exons(transcript)
		if len(exons) > 0:
			if exons[0].strand == "+":
				return exons[0]
			else:
				return exons[1]


	def get_transcript(self, exon):
				
		for transcript in self.exons[exon.chr][exon].keys():

			return transcript

	def get_transcripts(self, exon):
		return self.exons[exon.chr][exon].keys()

	def get_transcript_clusters(self):
		exons = self.get_exons()
		clusters = {}
		clus = {}
		c = 0
		for exon in exons:
			cluster = None
			for transcript in exon.evidence:
				if clusters.has_key(transcript):
					cluster = clusters[transcript]
			if not cluster:
				c += 1
				cluster = c
			for transcript in exon.evidence:
				clusters[transcript] = cluster
		for transcript,cluster in clusters.items():
			clus.setdefault(cluster, []).append(transcript)
		return clus.values()

	def get_chrs(self):
		return self.exons.keys()

	def get_leftmost_exons(self):
		return [exon for exon in self.get_exons() if not exon.linked]
	
	def get_valid_gene_models(self, exon):
		pass

	def get_models(self, exon):
		models = []
		for ev in exon.evidence:
			models += exon.linked_chains(ev)
		return models
	
	def get_exons_in_interval(self, chr, start, end, transcript=None):
		if not self.exons.has_key(chr):
			return []
		ret_exons = []
		for exon in self.exons[chr].keys():
			if (exon.start <= end and exon.start >= start) or (exon.end <= end and exon.end >= start):
				if not transcript or transcript in self.get_transcripts(exon):
					ret_exons.append(exon)
		return ret_exons

	def get_transcripts_in_interval(self, chr, start, end):
		transcripts = {}
		exons = self.get_exons_in_interval(chr, start, end)
		for exon in exons:
			for transcript in self.get_transcripts(exon):
				transcripts[transcript] = 1
		return transcripts.keys()

	def update_with(self,c, first_only=False):
		upstream = {}
		downstream = {}
		extra_exons = {}
		all_clusters =  c.get_transcript_clusters()
		for p, clus in enumerate(all_clusters):
			sys.stderr.write("%s of %s\n" % (p + 1, len(all_clusters)))
			if first_only:
				sorted_cluster = sorted(clus, key=lambda x: c.get_transcript_length(x))[:1]
			else:
				 sorted_cluster = sorted(clus, key=lambda x: c.get_transcript_length(x))
			for t in sorted_cluster:
				exons = c.get_transcript_exons(t)
				anno_exons = self.get_exons_in_interval(exons[0].chr, exons[0].start, exons[0].end)
				if anno_exons and anno_exons[-1].end == exons[0].end:
					exons[0].start = anno_exons[0].start
				
				anno_exons = self.get_exons_in_interval(exons[-1].chr, exons[-1].start, exons[-1].end)
				if anno_exons and anno_exons[0].start == exons[-1].start:
					exons[-1].end = anno_exons[0].end

				all_check = {}	
				for transcript in self.get_transcripts_in_interval(exons[0].chr, exons[0].start, exons[-1].end):	
					#print "Checking %s" % transcript	
					check = {}
					for exon in exons:
						update_exons = self.get_exons_in_interval(exon.chr, exon.start, exon.end, transcript=transcript)
						#print exon, update_exons
						if update_exons:
							if update_exons[0].start == exon.start and update_exons[-1].end == exon.end:
								check[exon] = 1
								all_check[exon] = 1
							else:
								check[exon] = 0
						else:
							check[exon] = 0


					if 1 in check.values():
						vals = [check[e] for e in sorted(check.keys(), key=lambda x: x.start)]
						#print vals
						if vals[0] == 0:
							i = 0
							while vals[i] == 0:
								i += 1
							t_exons = self.get_transcript_exons(transcript)
							
							if exons[0].end < t_exons[0].start and exons[0].strand == t_exons[0].strand: 
								
								if self.get_exon(exons[i].chr, exons[i].start, exons[i].end) in t_exons:
									j = t_exons.index(self.get_exon(exons[i].chr, exons[i].start, exons[i].end))
									self.delete_left_exons(transcript, j)
						
									for j in range(i):
										self.add(exons[j].chr, exons[j].start, exons[j].end, exons[j].strand, transcript=transcript)
	
									self.link_transcript_exons(transcript)
									#print self.get_transcript_exons(transcript)
									#print
									if exons[0].strand == "+":
										upstream[transcript] = 1
									else:
										downstream[transcript] = 1
							
						if vals[-1] == 0:
							i = len(vals) - 1
							while vals[i] == 0:
								i -= 1
							
							t_exons = self.get_transcript_exons(transcript)
							if exons[-1].start > t_exons[-1].end and exons[-1].strand == t_exons[-1].strand: 
									#print "exon in ",t_exons, self.get_exon(exons[i].chr, exons[i].start, exons[i].end)
								if self.get_exon(exons[i].chr, exons[i].start, exons[i].end) in t_exons:
									j = t_exons.index(self.get_exon(exons[i].chr, exons[i].start, exons[i].end))
									if len(t_exons) - j - 1 > 0:
										if exons[0].strand == "+":
											downstream[transcript] = 1
										else:
											upstream[transcript] = 1
										
										self.delete_right_exons(transcript, len(t_exons) - j - 1)
									
										for j in range(i + 1, len(vals)):
											self.add(exons[j].chr, exons[j].start, exons[j].end, exons[j].strand, transcript=transcript)
		
										self.link_transcript_exons(transcript)
										#print self.get_transcript_exons(transcript)
									#print
	
						for i in range(len(vals) - 2):
							if vals[i:i+3] == [1,0,1]:
								
								sys.stderr.write("Extra exon!\n")
								extra_exon = exons[i + 1]
								left = self.get_exons_in_interval(exons[i].chr, exons[i].start, exons[i].end, transcript)[-1]
								right = self.get_exons_in_interval(exons[i+2].chr, exons[i+2].start, exons[i+2].end, transcript)
								if right:
									right = right[0]
									between = self.get_exons_in_interval(left.chr, left.end + 1, right.start - 1, transcript)
								if between:
									sys.stderr.write("Complicated: %s\n" % extra_exon)
								else:
									#print "Easy"
							
									if left.links.has_key(right):
										transcripts = left.links[right]
										if transcript in transcripts:
											extra_exons[transcript] = 1
											#print "Before: ", self.get_transcript_exons(t)
											self.add(extra_exon.chr, extra_exon.start, extra_exon.end, extra_exon.strand, transcript=transcript)
											left.links[right].remove(transcript)
											if not left.links[right]:
												del left.links[right]
											self.link_exons(left, extra_exon, transcript)
											self.link_exons(extra_exon, right, transcript)

											#print "After: ", self.get_transcript_exons(t)


				if not 1 in all_check.values():
					if len(exons) > 1:
						for exon in exons:
							self.add(exon.chr, exon.start, exon.end, exon.strand, transcript="NEW_%s" % t)
					#print "all_new", exons
					
		return upstream, downstream, extra_exons


	def to_bed(self, model, name=None):	
		exons = self.get_transcript_exons(model)
		chromStart = exons[0].start
		chromEnd = exons[-1].end
		if not name:
			name = model
		sizes = ",".join([str(exon.end - exon.start) for exon in exons]) + ","
		starts = ",".join([str(exon.start - chromStart) for exon in exons]) + ","
		
		return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (exons[0].chr, chromStart, chromEnd, name, 600, exons[0].strand, chromStart, chromEnd, 0, len(exons), sizes, starts)


	def get_transcript_coords(self, transcript):
		exons = self.get_transcript_exons(transcript)
		#print "LALA",transcript,exons
		return [exons[0].chr, exons[0].start, exons[-1].end]

	def get_transcript_length(self, transcript):
		exons = self.get_transcript_exons(transcript)
		#print "LALA",transcript,exons
		return exons[-1].end -  exons[0].start

	def get_overlapping_transcripts(self, transcripts):
		#print len(transcripts)
		overlap = []
		ts = []
		for transcript in transcripts:
			ts.append([transcript] + self.get_transcript_coords(transcript))
		
		transcripts = sorted(ts, key=lambda x: (x[1], x[2]))
		i = 0
		while i < len(transcripts):
			t = [transcripts[i][0]]
			coords = transcripts[i][1:]
			#print coords
			i += 1
			if i < len(transcripts):
				coords2 = transcripts[i][1:]
				#print coords2
				while i < len(transcripts) and coords[0] == coords2[0] and coords2[1] >= coords[1] and coords2[1] <= coords[2]:
					t.append(transcripts[i][0])
					i += 1
					if i < len(transcripts):
						coords2 = transcripts[i][1:]
			overlap.append(t)
		#print "OVERLAP", overlap
		return overlap

	### Overlaps with external files (BAM, bed, etc.) ###
	def get_read_statistics(self, fname, name, span="exon", extend=(0,0)):
		from fluff.fluffio import get_binned_stats
		from tempfile import NamedTemporaryFile

		tmp = NamedTemporaryFile()
		estore = {}
		for exon in self.get_exons():
			start = exon.start
			end = exon.end
			if exon.strand == "-":
				start -= extend[1]
				end += extend[0]
			else:
				start -= extend[0]
				end += extend[1]
			if start < 0:
				start = 0
			
			estore["%s:%s-%s" % (exon.chr, start, end)] = exon
			tmp.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
				exon.chr,
				start,
				end,
				str(exon),
				0,
				exon.strand
			))
		tmp.flush()

		if fname.endswith("bam"):
			rmrepeats = True
		else:
			rmrepeats = False
		result = get_binned_stats(tmp.name, fname, 1, rpkm=False, rmdup=True, rmrepeats=rmrepeats)
		for row in result:
			vals = row.strip().split("\t")
			e = "%s:%s-%s" % (vals[0], vals[1], vals[2])
			c = float(vals[3])
			estore[e].stats[name] = c

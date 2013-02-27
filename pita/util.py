from solexatools import peak_stats
from solexatools.track import SimpleTrack
from tempfile import NamedTemporaryFile
import sys

def get_expression(exons):
	c = 0
	for exon in exons:
		print "EXON %s\t%s " %  (exon, exon.num_rnaseq)
		if exon.num_rnaseq >= MINIMUM_RNASEQ_READS:
			c += 1

	return float(c) / len(exons)

def same_exons(models):
	if len(models) == 1:
		return True
	
	exons = {}
	for model in models:
		exons[model] = sorted(c.get_transcript_exons(model), key=lambda x: x.start)

	l = len(exons[models[0]])
	for model in models[1:]:
		if len(exons[model]) != l:
			return False
	
	for i in range(len(exons[models[0]])):
		start, end = exons[models[0]][i].start, exons[models[0]][i].end
		for model in models[1:]:
			if exons[model][i].start != start and exons[model][i].end != end:
				return False
	return True
	
def get_closest_upstream(exon, pos):
	""" exon = Exon object, pos = sorted list (start, end) coordinates """
	if exon.strand == "+":
		i = 0
		while i < len(pos) and pos[i][0] < exon.end:
			i += 1
		if i > 0 and i < len(pos):
			return pos[i - 1]
	elif exon.strand == "-":
		i = len(pos) - 1
		while i >= 0 and pos[i][1] > exon.start:
			i -= 1
		if i + 1 < len(pos):
			if pos[i + 1][1] > exon.start:
				return pos[i + 1]
	else:
		raise ValueError, "Exon has no strand!"

def check_polII_upstream(d, polII_data):
	tmp = NamedTemporaryFile(delete=False)
	sys.stderr.write("check_polII_upstream: %s\n" % (tmp.name))
	for model, upstream in sorted(d.items()):
		model_coords = c.get_transcript_coords(model)
		tmp.write("%s\t%s\t%s\n" % tuple(model_coords))
		tmp.write("%s\t%s\t%s\n" % tuple(upstream))
	
	tmp.flush()
	result = peak_stats.peak_stats(SimpleTrack(tmp.name), SimpleTrack(polII_data), peak_stats.tuple_number_formatter, zeroes=True)
	true_upstream = []
	num = {}
	for row in result:
		num[(row[0], int(row[1]), int(row[2]))] = float(row[3])
	
	for (model, upstream) in sorted(d.items()):
		model_coords = c.get_transcript_coords(model)
		#sys.stderr.write("Checking model %s, %s versus upstream %s, %s vs %s\n" % (str(model), str(model_coords), str(upstream), num[tuple(model_coords)], num[tuple(upstream)] ))
		model_length = float(model_coords[2] - model_coords[1])
		up_length = float(upstream[2] - upstream[1])
		if (num[tuple(upstream)]/ up_length) / (( num[tuple(model_coords)] + 1)/model_length) > POLII_FRACTION_UPSTREAM:
			true_upstream.append(model)
	return true_upstream
			
def check_polII_at_tss(transcripts):
	
	dat = {}
	for transcript in transcripts:
		exons = c.get_transcript_exons(transcript)
	
		chr = exons[0].chr

		if exons[0].strand =="+":
			start = exons[0].start - 750
			end = exons[0].start + 250
		else:
			start = exons[-1].end - 250
			end = exons[-1].end + 750

		dat.setdefault((chr, start - 1000, start), []).append([transcript, "left"])
		dat.setdefault((chr, start, end), []).append([transcript, "middle"])
		dat.setdefault((chr, end, end + 1000), []).append([transcript, "right"])
		
	tmp = NamedTemporaryFile()
	for (chr, s, e) in dat.keys(): 
		tmp.write("%s\t%s\t%s\n" % (chr, s, e))
	
	tmp.flush()
	result = peak_stats.peak_stats(SimpleTrack(tmp.name), SimpleTrack(polII_data), peak_stats.tuple_number_formatter, zeroes=True)
	tmp.close()
	stats = {}
	for (chr, start, end, num) in result:
		for t,x in dat[(chr, start, end)]:
			stats.setdefault(t, {})[x] = float(num)
	
	ret = []
	for transcript in transcripts:
		if stats[transcript]["middle"] > stats[transcript]["left"] and stats[transcript]["middle"] > stats[transcript]["right"]:
			#print transcript, stats[transcript]["middle"], stats[transcript]["left"],  stats[transcript]["right"]
			ret.append(transcript)
	return ret
	
def get_updated_model(models):
	exons = {}
	all_exons = {}
	for model in models:
		exons[model] = sorted(c.get_transcript_exons(model), key=lambda x: x.start)
		for exon in exons[model]:
			all_exons[exon] = 1
	
	for exon in all_exons:
		if sorted(exon.evidence) != sorted(models):
			return True
	return False

def get_keys_of_highest_values(d, filter=None):
	if not filter:
		filter = d.keys()
	max = sorted([d[x] for x in d.keys() if x in filter])[-1]
	return [x for x in d.keys() if (d[x] == max) and ((not filter) or x in filter)]

def get_keys_of_lowest_values(d, filter=None):
	if not filter:
		filter = d.keys()
	min = sorted([d[x] for x in d.keys() if x in filter])[0]
	return [x for x in d.keys() if (d[x] == min) and ((not filter) or x in filter)]

def get_all_likely_models(models, status_id):
	if len(models) <= 1:
		return models

	overlap = []
	for i, model1 in enumerate(models):
		for model2 in models[i + 1:]:
			sys.stderr.write("%s\t%s\t%s\n" % ( model1, model2, get_exon_overlap(model1, model2)))
			o1,o2 = get_exon_overlap(model1, model2)
			overlap.append([o1, model1, model2])
			overlap.append([o2, model2, model1])

	filt_models = models[:]
	OVERLAP_THRESHOLD=0.8
	while len(overlap) > 1 and [x for x in overlap if x[0] >= OVERLAP_THRESHOLD]:
		key = sorted(overlap, key=lambda x: x[0])[-1][1]
		for i in range(len(overlap) - 1, -1, -1):
			if key in overlap[i]:
				del overlap[i]
		filt_models.remove(key)

	if len(filt_models) == 1:
		return filt_models
	
	tmp_models = filt_models[:]
	for model in filt_models[:]:
		if status_id[model].endswith("NX") or status_id[model].endswith("NU"):
			filt_models.remove(model)
	
	if len(filt_models) == 0:
		return [sorted(tmp_models, cmp=lambda x,y: cmp(len(c.get_transcript_exons(x)), len(c.get_transcript_exons(y))) )[-1]]


	# If we have removed everything, then they're all as likely
	#if len(filt_models) == 0:
	#	return models
	return filt_models

def get_most_likely_model(models, expression=True):
	""" Get the most likely model, based on expression, length, number of exons, etc"""
	# First, easy case: only one model
	if len(models) == 1:
		return models[0]

	if len(models) == 0:
		raise ValueError, "Length of models should not be 0"

	k4_models = [model for model in models if first_exon_has_k4(model)]
	if len(k4_models) == 1:
		return k4_models[0]
	elif len(k4_models) > 1:
		models = k4_models
	
	expressed_exons = {}
	number_of_exons = {}
	length = {}
	
	for model in models:
		exons = c.get_transcript_exons(model)
		number_of_exons[model] = len(exons)
		if expression:
			expressed_exons[model] = get_expression(exons) * len(exons)

		length[model] = exons[-1].end - exons[0].start
	#print expressed_exons
	if expression:
		likely_models = get_keys_of_highest_values(expressed_exons)
		if len(likely_models) == 1:
			return likely_models[0]
	else:
		likely_models = models[:]

	likely_models = get_keys_of_lowest_values(number_of_exons, filter=likely_models)
	if len(likely_models) == 1:
		return likely_models[0]
	
	if get_fm(likely_models):
		return get_fm(likely_models)
	
	likely_models = get_keys_of_highest_values(length, filter=likely_models)
	
	return likely_models[0]

def get_exon_overlap(t1, t2):
	exons1 = c.get_transcript_exons(t1)
	exons2 = c.get_transcript_exons(t2)
	
	overlap = 0.0
	l1 = 0
	for exon1 in exons1:
		l1 += len(exon1)
	l2 = 0
	for exon2 in exons2:
		l2 += len(exon2)
	
	for exon1 in exons1:
		for exon2 in exons2:
			if exon2.start >= exon1.start and exon2.start <= exon1.end:
				pos = exon1.end if exon2.end > exon1.end else exon2.end
				overlap += pos - exon2.start
			elif exon2.end >= exon1.start and exon2.end <= exon1.end:
				pos = exon1.start if exon1.start > exon2.start else exon2.start
				overlap += exon2.end - pos
	return overlap / l1, overlap / l2


def get_alternative_tss(models):
	#print "Check upstream TSS of %s" % str(models)
	for model in models:
		#if model == "ENSXETT00000034743":
		#	print 
		#	print "HOEI"
		#	print
#		print "LALA model: %s" % model
		exons = c.get_transcript_exons(model)
		extra = None
		gene = [exons[0].chr, exons[0].start, exons[-1].end]
		if exons[0].strand == "-" and pos.has_key(exons[-1].chr):
			p =  get_closest_upstream(exons[-1], pos[exons[-1].chr])
#			print "LALA peak: %s" % str(p)
			if p:
				p_start, p_end = p
				same_loc_exons = sorted(c.get_exons_in_interval(exons[-1].chr, exons[-1].start, exons[-1].end), key=lambda x: x.start)
				exons_in_between = c.get_exons_in_interval(exons[-1].chr, exons[-1].end, p_start)
				transcripts_in_between = {}
				for exon in exons_in_between:
					if exon.start > exons[-1].end:
						for t in c.get_transcripts(exon):
							if t not in models and not t.startswith("NEW"):
								transcripts_in_between[t] = 1
				
#				print "LALA in between: %s" % str(transcripts_in_between.keys())
				
				if len(transcripts_in_between.keys()) == 0:
					extra = [exons[0].chr, exons[-1].end, p_start]

		elif exons[0].strand == "+" and pos.has_key(exons[-1].chr):
			p =  get_closest_upstream(exons[0], pos[exons[-1].chr])
			if p:
				p_start, p_end = p
				same_loc_exons = sorted(c.get_exons_in_interval(exons[0].chr, exons[0].start, exons[0].end), key=lambda x: x.start)
				lala = same_loc_exons[0].start - 1	
				exons_in_between = c.get_exons_in_interval(exons[-1].chr, p_end, same_loc_exons[0].start - 1)
				transcripts_in_between = {}
				for exon in exons_in_between:
					if exon.end < exons[0].start:
						for t in c.get_transcripts(exon):
							if t not in models and not t.startswith("NEW"):
								transcripts_in_between[t] = 1

				if len(transcripts_in_between.keys()) == 0:
					extra = [exons[0].chr, p_end, exons[0].start]
		#print extra
		#print
		if extra:
			#return "Should be checked, xtra: %s:%s-%s gene: %s:%s-%s" % (extra[0], extra[1], extra[2], gene[0], gene[1], gene[2])
			return TSS_UPSTREAM, (extra[0], extra[1], extra[2])
		
	return TSS_NOTFOUND, None


def get_downstream_tss(transcript):
	tss = []
	s_exons = c.get_transcript_exons(transcript)

	if s_exons[0].strand == "+":
		tss_k4 = s_exons[0].k4
		exons = s_exons[1:]
	else:
		tss_k4 = s_exons[-1].k4
		exons = s_exons[:-1]
	
	for exon in exons:
		if exon.k4 and exon.k4[0] not in tss_k4:
			if len(tss_k4) == 0 or float(exon.k4[0].split("\t")[3]) > 0.3 * float(tss_k4[0].split("\t")[3]):
				#print "LALA", float(exon.k4[0].split("\t")[3]),  0.3 * float(tss_k4[0].split("\t")[3])
				tss_k4.append(exon.k4[0])
				tss.append(exon)
	
	return tss


def get_expressed_models(models):
	print "GET expressed models"
	for t in models:
		print "%s\t%s" % (t,  get_expression(c.get_transcript_exons(t)))
	return [transcript for transcript in models if get_expression(c.get_transcript_exons(transcript)) >= MINIMUM_FRACTION_EXPRESSED]

def first_exon_has_k4(transcript):
	s_exons = c.get_transcript_exons(transcript)
	if s_exons[0].strand == "+":
		#print "check k4 %s\t%s" % (transcript, str(s_exons[0].k4))
		return s_exons[0].k4
	else:
		#print "check k4 %s\t%s" % (transcript, str(s_exons[-1].k4))
		return s_exons[-1].k4

def get_tss_models(models):
	return [transcript for transcript in models if first_exon_has_k4(transcript)]

#def get_fm(models):
#	for model in models:
#		try:
#			int(model)
#			return model
#		except:
#			pass
def get_fm(models):
	for model in models:
		if model.startswith("xb_"):
			return model

def get_k4_peaks(tmpname, k4_peaks):
	sys.stderr.write("Determining H3K4me3 peak overlap\n")
	exon_lines = []
	s =  SimpleTrack(tmpname)
	f = s.get_next_feature()
	while f:
		exon_lines.append(f[:3])
		f = s.get_next_feature()
		
	sys.stderr.write("Determining H3K4me3 peak overlap (1)\n")
	result = peak_stats.peak_stats(SimpleTrack(tmpname), SimpleTrack(k4_peaks), peak_stats.tuple_all_formatter, zeroes=True)
	sys.stderr.write("Determining H3K4me3 peak overlap (2)\n")
	
	for row,exon in zip(result, exon_lines):
		chr,start,end = exon
		start,end = int(start),int(end)
		for strand in ["+", "-"]:
			exon = c.get_exon(chr, start, end, strand)
			if exon:
				exon.k4 = row[:]	
	sys.stderr.write("Determining H3K4me3 peak overlap (done)\n")

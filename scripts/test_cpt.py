#!/usr/bin/env python
from pita.io import *
from pita.util import bed2exonbed
from tempfile import NamedTemporaryFile
import subprocess as sp
import sys
import os
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import numpy as np
import argparse

cp = importr('changepoint')

def call_cpt(start, end, strand, data, min_reads=5, min_log2_ratio=1.5):
    pt_cutoff = 5
    
    #sys.stderr.write("{}\t{}\t{}\t{}\n".format(start, end, strand, len(data)))
    counts = np.array(data)
    #print counts 
    if strand == "-":
        counts = counts[::-1]
    r_counts = robjects.FloatVector([float(c) for c in counts])
    result = cp.cpt_mean(r_counts)
    try:
        pt = int(cp.cpts(result)[0])
    except:
        pt = 0
        
    #sys.stderr.write("{}-{}: cpt is {}\n".format(start, end, pt))
    if pt > pt_cutoff:
        ratio = np.log2(counts[:pt].mean() / (counts[pt:].mean()) + 0.1)
    while pt > pt_cutoff and ratio < 1:
        r_counts = robjects.FloatVector(counts[:pt])
        result = cp.cpt_mean(r_counts)
        try: 
            pt = int(cp.cpts(result)[0])
        except:
            pt = 0
        #sys.stderr.write("{}-{}: updating cpt to {}\n".format(start, end, pt))
        if pt > pt_cutoff:
            ratio = np.log2(counts[:pt].mean() / (counts[pt:].mean()) + 0.1)
        
    if pt > pt_cutoff:
        while pt < len(counts) and counts[pt] >= min_reads:
            pt += 1
            
        m = counts[:pt].mean()
        md = np.median(counts[:pt])
        s = np.std(counts[:pt])
        q = np.percentile(counts[:pt], 0.25)
        
        if m >= min_reads and ratio >= min_log2_ratio:
            
            if strand == "-":
                utr_start = int(end) - pt
                utr_end = int(end)
            else:
                utr_start = int(start)
                utr_end = int(start) + pt
    
            return utr_start, utr_end
           
def call_utr(inbed, bamfiles):
    sys.stderr.write("Preparing temporary BED files\n")
    exonbed = NamedTemporaryFile()
    bed2exonbed(inbed, exonbed.name)

    transcripts = read_bed_transcripts(open(bedfile))
    td = dict([(t[0].split("_")[1] + "_", t[2]) for t in transcripts])
    
    genes = pybedtools.BedTool(bedfile)
    exons = pybedtools.BedTool(exonbed.name)
    
    tmp = NamedTemporaryFile()
    
    EXTEND = 10000
    sys.stderr.write("Determining gene boundaries determined by closest gene\n")
    for x in genes.closest(exons, D="a", io=True, iu=True):
        transcript = td[x[3]]
        extend = EXTEND
        if (int(x[-1]) >= 0) and (int(x[-1]) < extend):
            extend = int(x[-1])
        
        if transcript[0][-1] == "+":
            first = transcript[-1]
            first[2] += extend
        else:
            first = transcript[-0]
            first[1] -= extend
    
            if first[1] < 0:
                first[1] = 0
       
        tmp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            first[0],
            first[1],
            first[2],
            x[3],
            0,
            first[3]
            ))
        
    tmp.flush()
    
    tmpsam = NamedTemporaryFile()
    tmpbam = NamedTemporaryFile()
    sp.call("samtools view -H {} > {}".format(bamfiles[0], tmpsam.name), shell=True)
    cmd = "samtools view -L {} {} >> {}"
    sys.stderr.write("Merging bam files\n")
    for bamfile in bamfiles:
        sp.call(cmd.format(tmp.name, bamfile, tmpsam.name), shell=True)
    cmd = "samtools view -Sb {} | samtools sort - {}"
    sp.call(cmd.format(tmpsam.name, tmpbam.name), shell=True)
    sp.call("samtools index {}.bam".format(tmpbam.name), shell=True)
    tmpsam.close()
    
    sys.stderr.write("Calculating coverage\n")
    cmd = "bedtools coverage -abam {} -b {} -d -split "
    
    p = sp.Popen(cmd.format(tmpbam.name + ".bam", tmp.name), shell=True, stdout=sp.PIPE, bufsize=1)
    
    sys.stderr.write("Calling UTRs\n")
    
    data = []
    current = [None]
    utr = {}
    for line in iter(p.stdout.readline, b''):
        vals = line.strip().split("\t")
        if vals[3] != current[0]:
            if len(data) > 0:
                result = call_cpt(current[1], current[2], current[3], data, len(bamfile))
                #print result
                if result:
                    utr[current[0]] = result
            data = []
            current = [vals[3], int(vals[1]), int(vals[2]), vals[5]]
        data.append(int(vals[7]))
    result = call_cpt(current[1], current[2], current[3], data, len(bamfiles))
    #print result
    if result:
        utr[current[0]] = result
    
    tmpbam.close()
    return utr

def print_updated_bed(bedfile, bamfiles):
    utr = call_utr(bedfile, bamfiles)
    for line in open(bedfile):
        
        vals = line.strip().split("\t")
        start,end = int(vals[1]), int(vals[2])
        strand = vals[5]
        name = vals[3]
        thickstart, thickend = int(vals[6]), int(vals[7])
        exonsizes = [int(x) for x in vals[10].split(",") if x]
        exonstarts = [int(x) for x in vals[11].split(",") if x]
        
        if utr.has_key(name):
            sys.stderr.write("Updating {}\n".format(name))
            utr_start, utr_end = utr[name]
        
            if strand == "+":
                if utr_end < thickend:
                    sys.stderr.write("Setting end of {} to CDS end\n".format(name))
                    utr_end = thickend
                diff = exonsizes[-1] - (utr_end - utr_start)
                end -= diff
                
                exonsizes[-1] -= diff
    
                vals[2] = end
                vals[10] = ",".join([str(x) for x in exonsizes] + [""])
            else:
                if utr_start > thickstart:
                    sys.stderr.write("Setting start of {} to CDS start\n".format(name))
                    utr_start = thickstart
                diff = exonsizes[0] - (utr_end - utr_start)
                sys.stderr.write("{} {} {} diff: {}\n".format(utr_start, utr_end, exonsizes[0], diff))
                start += diff
                
                exonstarts = [0] + [x - diff for x in exonstarts[1:]]
                exonsizes[0] -= diff
    
                vals[1] = start
                vals[10] = ",".join([str(x) for x in exonsizes] + [""])
                vals[11] = ",".join([str(x) for x in exonstarts] + [""])
    
            print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(*vals)
        else:
            print line.strip()

p = argparse.ArgumentParser()
p.add_argument("-i",
               dest= "bedfile",
               help="genes in BED12 format",
              )
p.add_argument("-b",
               dest= "bamfiles",
               help="list of RNA-seq BAM files (seperated by comma)",
              )
args = p.parse_args()

if not args.bedfile or not args.bamfiles:
    p.print_help()
    sys.exit()

bedfile = args.bedfile
bamfiles = args.bamfiles.split(",")

for bamfile in bamfiles:
    if not os.path.exists(bamfile):
        sys.stderr.write("BAM file {} does not exist.\n".format(bamfile))
        sys.exit(1)    
    if not os.path.exists(bamfile + ".bai"):
        sys.stderr.write("index file {}.bai does not exist.\n".format(bamfile))
        sys.exit(1)

if not os.path.exists(bedfile):
    sys.stderr.write("BED file {} does not exist.\n".format(bedfile))
    sys.exit(1)

print_updated_bed(bedfile, bamfiles)

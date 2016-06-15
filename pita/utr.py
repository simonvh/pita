from pita.r_cpt import cpt
from pita.util import bed2exonbed
from pita.io import read_bed_transcripts
from tempfile import NamedTemporaryFile
import pybedtools
import subprocess as sp
import sys
import os
import numpy as np

def call_cpt(start, end, strand, data, min_reads=5, min_log2_ratio=1.5, upstream=False):
    """
    Determine UTR location from basepair-resolution coverage vector of reads.
    Return tuple of utr_start and utr_end if changepoint is found.
    """
    
    sys.stderr.write("{} {} {}\n".format(start, end, strand))
    pt_cutoff = 5
    counts = np.array(data)
    
    # Do calculations on reverse array if gene is on the - strand or
    # when predicting 5' UTR
    if (upstream and strand == "+") or strand == "-":
        counts = counts[::-1]
    
    pt = len(counts)
    ratio = 0
    while pt > pt_cutoff and ratio < 1:
        sys.stderr.write("cpt {}\n".format(pt))
        pt = cpt(counts[:pt])
        if pt > pt_cutoff:
            ratio = np.log2(counts[:pt].mean() / (counts[pt:].mean()) + 0.1)
        
    if pt > pt_cutoff:
        # Add to the changepoint while the number of reads is above min_reads
        while pt < len(counts) and counts[pt] >= min_reads:
            pt += 1
            
        m = counts[:pt].mean()
        if m >= min_reads and ratio >= min_log2_ratio:
            if (upstream and strand == "+") or strand == "-":
                utr_start = int(end) - pt
                utr_end = int(end)
            else:
                utr_start = int(start)
                utr_end = int(start) + pt
    
            return int(utr_start), int(utr_end)
           
def call_utr(inbed, bamfiles, utr5=False, utr3=True):
    """
    Call 3' UTR for all genes in a BED12 file based on RNA-seq reads 
    in BAM files.
    """
    
    # Load genes in BED file
    transcripts = read_bed_transcripts(open(inbed))
    
    # No genes
    if len(transcripts) == 0:
        return 

    td = dict([(t[0].split("|")[1] + "_", t[2]) for t in transcripts])
    
    #Trying to fix the scaffold struggles
    #td = {}
    #for t in transcripts:
    #    if "scaffold" not in t[0]:
    #        td[t[0].split("_")[1]+"_"] = t[2]
#	else:
#	    inter = t[0].split(":")
#            scafName = "_".join(inter[0].split("_")[2:4])
#	    pos = inter[1].split("_")[0]
#	    td[scafName+":"+pos+"_"] = t[2]

    # Create a BED6 file with exons, used to determine UTR boundaries 
    sys.stderr.write("Preparing temporary BED files\n")
    exonbed = NamedTemporaryFile(prefix="pita.", suffix=".bed")
    bed2exonbed(inbed, exonbed.name)

    # Determine boundaries using bedtools
    genes = pybedtools.BedTool(inbed)
    exons = pybedtools.BedTool(exonbed.name)
    
    tmp = NamedTemporaryFile(prefix="pita.", suffix=".bed")
    
    EXTEND = 10000
    sys.stderr.write("Determining gene boundaries determined by closest gene\n")
    for x in genes.closest(exons, D="a", io=True, iu=True):
        transcript = td[x[3]]
        
        # Extend to closest exon or EXTEND, whichever is closer
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
    
    tmpsam = NamedTemporaryFile(prefix="pita.", suffix=".sam")
    tmpbam = NamedTemporaryFile(prefix="pita.")
    
    # Retrieve header from first BAM file
    sp.call("samtools view -H {} > {}".format(bamfiles[0], tmpsam.name), shell=True)
    
    # Filter all BAM files for the specific regions. This runs much faster
    # then running bedtools coverage on all individual BAM files
    tmp_check = NamedTemporaryFile(prefix="pita.", suffix=".bam")
    cmd = "samtools view -L {} {} > {}"
    sys.stderr.write("Merging bam files\n")
    for bamfile in bamfiles:
        try:
            sp.check_call(cmd.format(tmp.name, bamfile, tmp_check.name), shell=True)
            sp.call("cat {} >> {}".format(tmp_check.name, tmpsam.name), shell=True)
        except sp.CalledProcessError as e:
            sys.stderr.write("Error in file {}, skipping:\n".format(bamfile))
            sys.stderr.write("{}\n".format(e))
    
    tmp_check.close()

    # Created sorted and index bam
<<<<<<< HEAD
    cmd = "samtools view -Sb {} | samtools sort -m 4G - {}"
=======
    cmd = "samtools view -Sb {} | samtools sort -m 10G - {}"
>>>>>>> master
    sp.call(cmd.format(tmpsam.name, tmpbam.name), shell=True)
    sp.call("samtools index {}.bam".format(tmpbam.name), shell=True)
    
    # Close and remove temporary SAM file
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
    if current[0]:
        result = call_cpt(current[1], current[2], current[3], data, len(bamfiles))
        if result:
            utr[current[0]] = result
    
    for fname in [tmpbam.name + ".bam", tmpbam.name + ".bam.bai"]:
        if os.path.exists(fname):
            os.unlink(fname)
    
    tmpbam.close()
    tmp.close()
    
    
    return utr

# A "nice" hack to implement 5' en 3' utr extension
def flip_bed_strands(bedfile):
    temp = NamedTemporaryFile(delete=False)
    for line in open(bedfile):
        line = line.strip().split("\t")
        if line[5] == "-":
            line[5] = "+"
        elif line[5] == "+":
            line[5] = "-"
        temp.write("\t".join(line)+"\n")
    temp.flush()
    return temp.name

def print_updated_bed(bedfile, bamfiles):
    #Extend the utr to the 5' end
    first = calculate_updated_bed(bedfile, bamfiles)

    #flipping the strands to extend the utr to the 3' end
    preSecond = flip_bed_strands(first)
    second = calculate_updated_bed(preSecond, bamfiles)

    #back to the correct strands
    final = flip_bed_strands(second)
    
    #print the utrExtended bed to the console
    for line in open(final):
        print(line.strip())


def calculate_updated_bed(bedfile, bamfiles):
    temp = NamedTemporaryFile(delete=False)

    utr = call_utr(bedfile, bamfiles)
    for line in open(bedfile):
        if line.startswith("track") or line[0] == "#":
            temp.write(line.strip())
            continue
        
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
    
            temp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*vals))
        else:
            temp.write(line)
    temp.flush()
    return temp.name


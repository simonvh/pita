from pita.config import SAMTOOLS
from pita.config import PitaConfig
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
import sys
import subprocess as sp
from tempfile import NamedTemporaryFile
import logging
import os
import sys
import yaml
import argparse

logger = logging.getLogger('pita')

def read_statistics(fname, rmrepeat=False, rmdup=False, mapped=False):
    """ Count number of reads in BAM file.
    Optional arguments rmrepeat and rmdup do nothing for now
    """

    cmd = "{0} idxstats {1} | awk '{{total += $3 + $4}} END {{print total}}'"

    p = Popen(
             cmd.format(SAMTOOLS, fname), 
             shell=True, 
             stdout=PIPE, 
             stderr=PIPE
             )

    stdout,stderr = p.communicate()

    n = int(stdout.strip())

    return n

def get_overlapping_models(exons):
    overlap = []
    sorted_exons = sorted(exons, cmp=lambda x,y: cmp(x.start, y.start))

    for i, exon in enumerate(sorted_exons):
        j = i + 1
        while j < len(sorted_exons) and sorted_exons[j].start <= exon.end:
            if exon.overlap(sorted_exons[j]):
                overlap.append([exon, sorted_exons[j]])
            j += 1
    
    return overlap

def exons_to_seq(exons):
    seq = ""
    if exons[0].strand == "-":
        exons = exons[::-1]

    for e in exons:
        if e.seq:
            seq = "".join((seq, e.seq))
        else:
            logger.error("exon {} has no sequence".format(e))

    return seq

def longest_orf(seq, do_prot=False):
    if type(seq) == type([]):
        seq = exons_to_seq(seq)

    dna = Seq(seq, IUPAC.ambiguous_dna)

    my_cmp = lambda x,y: cmp(len(x), len(y))

    orfs = []
    prots = []
    for i in range(3):
        seq = dna[i:]

        # BioPython doesn't like translating DNA with a length that's not
        # a multiple of 3
        if len(seq) % 3:
            seq = seq[:-(len(seq) % 3)]

        prot = str(seq.translate())
        #putative_orfs = [re.sub(r'^[^M]*', "", o) for o in prot.split("*")]
        putative_orfs = [o for o in prot.split("*")]
        longest = sorted(putative_orfs, cmp=my_cmp)[-1]
        prots.append(longest)
        start = prot.find(longest) * 3 + i
        end = start + (len(longest) + 1) * 3
        orfs.append((start,end))

    if do_prot:
        return sorted(prots, cmp=lambda x,y: cmp(len(x), len(y)))[-1]
    else:
        return sorted(orfs, cmp=lambda x,y: cmp(x[1] - x[0], y[1] - y[0]))[-1]

def find_genomic_pos(pos, exons):
    if exons[0].strand == "-":
        for e in exons[::-1]:
            size = e.end - e.start
            if pos < size:
                return e.end - pos
            else:
                pos -= size
        return e.start
    else:
        for e in exons:
            size = e.end - e.start
            if pos < size:
                return e.start + pos
            else:
                pos -= size
        return e.end

def to_genomic_orf(start, end, exons):
    genomic_start = find_genomic_pos(start, exons)
    genomic_end = find_genomic_pos(end,  exons)
    if exons[0].strand == "-":
        return genomic_end, genomic_start
    else:
        return genomic_start, genomic_end

def model_to_bed(exons, genename=None):
    line_format = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}"
    
    chrom = exons[0].chrom
    chromStart = exons[0].start
    chromEnd = exons[-1].end
    
    if not genename:
        genename = "{0}:{1}-{2}".format(chrom, chromStart, chromEnd)
    
    thickStart, thickEnd = chromStart, chromEnd

    if exons[0].seq:
        orf_start, orf_end = longest_orf(exons)
        thickStart, thickEnd = to_genomic_orf(orf_start, orf_end, exons)

    sizes = ",".join([str(exon.end - exon.start) for exon in exons]) + ","
    starts = ",".join([str(exon.start - chromStart) for exon in exons]) + ","

    return line_format.format(
                            chrom,
                            chromStart,
                            chromEnd,
                            genename,
                            600,
                            exons[0].strand,
                            thickStart,
                            thickEnd,
                            "0,0,0",
                            len(exons),
                            sizes,
                            starts
                            )



def get_splice_score(a, s_type=5):
    if not s_type in [3,5]:
        raise Exception("Invalid splice type {}, should be 3 or 5".format(s_type))
    maxent = getMaxPath()
    tmp = NamedTemporaryFile()
    for name,seq in a:
        tmp.write(">{}\n{}\n".format(name,seq))
    tmp.flush()
    cmd = "perl score{}.pl {}".format(s_type, tmp.name)
    p = sp.Popen(cmd, shell=True, cwd=maxent, stdout=sp.PIPE)
    score = 0
    for line in p.stdout.readlines():
        vals = line.strip().split("\t")
        if len(vals) > 1:
            score += float(vals[-1])
    return score

#parse config file for the maxentPath
def getMaxPath():
	DEFAULT_CONFIG = "pita.yaml"
	DEFAULT_THREADS = 4
	p = argparse.ArgumentParser()
	p.add_argument("-c",
               	dest= "configfile",
               	default = DEFAULT_CONFIG,
               	help="configuration file (default:{0})".format(DEFAULT_CONFIG)
              )
	p.add_argument("-t",
               	dest= "threads",
               	default = DEFAULT_THREADS,
               	type = int,
               	help="number of threads (default:{0})".format(DEFAULT_THREADS)
              )
	p.add_argument("-i",
               	dest= "index_dir",
               	default = None,
               	help="genome index dir"
              )
	p.add_argument("-y",
               	dest= "yaml_file",
               	default = None,
               	help="dump database to yaml file"
              )
	p.add_argument("-d",
               	dest= "debug_level",
               	default = "INFO",
               	help="Debug level"
              )
	p.add_argument("-r",
 		dest ="reannotate",
                default = False,
                action = "store_true",
                help="reannotate using existing database"
                )
	args = p.parse_args()
	configfile = args.configfile
	if not os.path.exists(configfile):
		print("Configfile not found")
		sys.exit()
	config = PitaConfig(configfile,args.reannotate)

        return config.maxentpath

def bed2exonbed(inbed, outbed):
    with open(outbed, "w") as out:
        for line in open(inbed):
            if line.startswith("#") or line.startswith("track"):
                out.write(line)
            else:
                vals = line.strip().split("\t")
        
                exonsizes = [int(x) for x in vals[10].split(",") if x]
                exonstarts = [int(x) for x in vals[11].split(",") if x]
                for exon_start, exon_size in zip(exonstarts, exonsizes):
                    out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (vals[0], int(vals[1]) + exon_start, int(vals[1]) + exon_start + exon_size, vals[3], vals[4], vals[5]))
    

if __name__ == "__main__":
    print longest_orf("CAGGAAGTCACGGAGCGCGGGATTTTTCAATCAGACTGATGAACAGATGAATACGACGAAGAGCATGGAGGCAATTCTGGAATTTTTTGTGCTGTGTGATCCAAAGAAGCGGCCAGTCAGACTGAACCGGTTGCCTTCTGTACCAAAGGATGCACTGTGTTATTCTGCCCTGCTGCCATCTCCTCTACCATCCCAGCTGTTGATCTTTGGCTTAGGTGACTGGTCAGGGTTATCTGGAGGAAGCACAGTAGAAGTGAAATTGGAAGGAAGTGGAACCAAAGAGCACAGACTGGGAACGCTGACTCCTGAGTCAAGATGCTTCCTGTGGGAATCTGACCAAAACCCCGACACCAGCATAATGTTACAAGAGGGAAAGCTGCATATCTGCATGTCGGTTAAAGGGCAGGTCAATATTAATTCTACTAACAGGAAAAAAGAGCATGGAAAGCGCAAGAGAATTAAAGAGGAAGAGGAAAATGTTTGTCCAAATAGTGGACATGTAAAAGTGCCTGCTCAAAAACAGAAGAACAGTAGTCCTAAGAGTCCAGCACCAGCAAAGCAACTTGCTCATTCTAAGGCCTTTTTAGCAGCACCAGCTGTGCCAACTGCACGCTGGGGTCAAGCGCTCTGTCCTGTCAACTCTGAGACAGTAATCTTGATTGGTGGACAGGGAACACGTATGCAGTTCTGTAAGGATTCCATGTGGAAACTGAATACAGATAGGAGCACATGGACTCCAGCTGAGGCATTGGCAGATGGCCTTTCACCAGAAGCTCGTACTGGGCACACAGCAACCTTCGATCCTGAGAACAACCGTATTTATGTGTTTGGAGGTTCTAAGAACAGAAAATGGTTCAATGATGTACATATTTTGGACATTGAGGCCTGGCGATGGAGGAGCGTGGAAGTAAGTAAACTAAGTAGTTGA")

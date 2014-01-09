from pita.config import SAMTOOLS
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re

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
        if not e.seq:
            raise Exception, "exon has no sequence"
        seq = "".join((seq, e.seq))
    return seq


def longest_orf(seq):
    if type(seq) == type([]):
        seq = exons_to_seq(seq)
    dna = Seq(seq, IUPAC.ambiguous_dna)
    
    my_cmp = lambda x,y: cmp(len(x), len(y))

    orfs = []
    for i in range(3):
        prot = str(dna[i:].translate())
        putative_orfs = [re.sub(r'^[^M]*', "", o) for o in prot.split("*")]
        longest_orf = sorted(putative_orfs, cmp=my_cmp)[-1]
        start = prot.find(longest_orf) * 3 + i
        end = start + (len(longest_orf) + 1) * 3
        orfs.append((start,end))
    
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

if __name__ == "__main__":
    print longest_orf("CAGGAAGTCACGGAGCGCGGGATTTTTCAATCAGACTGATGAACAGATGAATACGACGAAGAGCATGGAGGCAATTCTGGAATTTTTTGTGCTGTGTGATCCAAAGAAGCGGCCAGTCAGACTGAACCGGTTGCCTTCTGTACCAAAGGATGCACTGTGTTATTCTGCCCTGCTGCCATCTCCTCTACCATCCCAGCTGTTGATCTTTGGCTTAGGTGACTGGTCAGGGTTATCTGGAGGAAGCACAGTAGAAGTGAAATTGGAAGGAAGTGGAACCAAAGAGCACAGACTGGGAACGCTGACTCCTGAGTCAAGATGCTTCCTGTGGGAATCTGACCAAAACCCCGACACCAGCATAATGTTACAAGAGGGAAAGCTGCATATCTGCATGTCGGTTAAAGGGCAGGTCAATATTAATTCTACTAACAGGAAAAAAGAGCATGGAAAGCGCAAGAGAATTAAAGAGGAAGAGGAAAATGTTTGTCCAAATAGTGGACATGTAAAAGTGCCTGCTCAAAAACAGAAGAACAGTAGTCCTAAGAGTCCAGCACCAGCAAAGCAACTTGCTCATTCTAAGGCCTTTTTAGCAGCACCAGCTGTGCCAACTGCACGCTGGGGTCAAGCGCTCTGTCCTGTCAACTCTGAGACAGTAATCTTGATTGGTGGACAGGGAACACGTATGCAGTTCTGTAAGGATTCCATGTGGAAACTGAATACAGATAGGAGCACATGGACTCCAGCTGAGGCATTGGCAGATGGCCTTTCACCAGAAGCTCGTACTGGGCACACAGCAACCTTCGATCCTGAGAACAACCGTATTTATGTGTTTGGAGGTTCTAAGAACAGAAAATGGTTCAATGATGTACATATTTTGGACATTGAGGCCTGGCGATGGAGGAGCGTGGAAGTAAGTAAACTAAGTAGTTGA")

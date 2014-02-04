import sys
import HTSeq
from pita.io import read_bed_transcripts
import numpy as np

def compare_annotation(models1, models2):

    ga = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    l = {}
    result = {}
    for name, x, exons in models1:
        #exons[0][1] = exons[0][2] - 10
        #exons[-1][2] = exons[-1][1] + 10
        l[name] = 0.0
        for exon in exons:
            iv = HTSeq.GenomicInterval(*exon)
            ga[iv] += name
            l[name] += iv.end - iv.start

    for name, x, exons in models2:
        #exons[0][1] = exons[0][2] - 10
        #exons[-1][2] = exons[-1][1] + 10
        overlap = {}
        tlen = 0.0
        for exon in exons:
            test = HTSeq.GenomicInterval(*exon)
            
            for iv2, step_set in ga[test].steps():
                if len(step_set) > 0:
                    for n in step_set:
                        overlap[n] = overlap.setdefault(n, 0) + iv2.end - iv2.start
            tlen += exon[2] - exon[1]
        
        acs = []
        
        result[name] = {}
        for oname, length in overlap.items():
       
            #print "#", name, tlen, oname, l[oname], length
            #print name, oname 
            sn = length / l[oname]
            sp = length / tlen
            ac = (sn + sp) / 2
            acs.append(ac)
            result[name][oname] = [sn, sp, ac]
    return result        
        
def compare_annotation_files(fname1, fname2):        
    file1 = open(fname1)
    file2 = open(fname2)

    data1 = read_bed_transcripts(file1)
    data2 = read_bed_transcripts(file2)
       
    result = compare_annotation(data1, data2)
    for t2, x, exons in data2:
        if not result.has_key(t2) or result[t2] == {}:
            print "{0}\t{1}\t{2}\t{3}".format(exons[0][0], exons[0][1], exons[-1][2], 1)

    #print result 
    cluster = {}
    t2_cluster = {}
    for transcript2,d in result.items():
        for transcript1 in d.keys():
            if not cluster.has_key(transcript1):
                cluster[transcript1] = d.keys()
            else:
                for t in d.keys(): 
                    if not t in cluster[transcript1]:
                        cluster[transcript1].append(t) 
            for t1 in cluster[transcript1]:
                t2_cluster.setdefault(t1, []).append(transcript2)
    
    clusters = set([tuple(x) for x in cluster.values()])
    
    loc = {}
    for name,x,exons in data1:
        loc[name] = [exons[0][0], exons[0][1], exons[-1][2]]
    
    for cluster in clusters:
        t2s = set(t2_cluster[cluster[0]])
        min_acs = []
        for transcript1 in cluster:
            acs = []
            for transcript2 in t2s:
                if result[transcript2].has_key(transcript1):
                    sn, sp, ac = result[transcript2][transcript1]
                    #print transcript1, transcript2, sn, sp, ac
                    acs.append(1 - sp)
            min_acs.append(np.min(acs))
            #print transcript1, np.min(acs)
        #print "###"
        #print cluster, min_acs
        chrom = loc[cluster[0]][0]
        start = np.min([loc[t][1] for t in cluster])
        end = np.max([loc[t][2] for t in cluster])
        print "{0}\t{1}\t{2}\t{3}".format(chrom, start, end, np.mean(min_acs))
    #print loc   

        
        #for t in cluster:
    sys.exit()
    
    for k,v in result.items():
        if len(v.keys()) > 0:
            print "{0}\t{1}\t{2}".format(loc[k], k, 1 - np.max(v.values()))
        else:
            print "{0}\t{1}\t{2}".format(loc[k], k, 1.0)
        











fname1 = "/home/simon/prj/laevis/annotation/XENLA_JGI7b/gmap/xlaevisEST.gmap.bed"  
#fname2 = "/home/simon/prj/laevis/annotation/XENLA_JGI7b/gmap/XENLA_2013may.longest_cdna_annot.XENLA_JGIv7b.gmap.bed"
fname2 = "/home/simon/prj/laevis/pita/v1.32/pita_v1.32.bed"


#x = HTSeq.GenomicInterval( "chr3", 123203, 127245, "+" )
#sys.exit()


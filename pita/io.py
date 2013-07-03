from BCBio import GFF
from pita.collection import *
import pprint
import sys
import logging 

# Setup logging
logger = logging.getLogger('pita')

def merge_exons(exons):
    return exons
    if not 0 in [e2.start - e1.end for e1,e2 in zip(exons[0:-1], exons[1:])]:
        return exons
    new_exons = []
    #for e2, e1

def _gff_type_iterator(feature, ftypes):
    if feature.type in ftypes:
        yield feature
    else:
        for feature in feature.sub_features:
            for f in _gff_type_iterator(feature, ftypes):
                yield f

def read_gff_transcripts(fname, min_exons=1):
    #limits = dict(gff_type = ["mRNA", "exon"])
    smap = {"1":"+",1:"+","-1":"-",-1:"-", None:"+"}
    transcripts = []
    for rec in GFF.parse(open(fname)):
        chrom = rec.id
        for feature in rec.features:
            logger.debug("feature: {0}".format(feature))
            
            for gene in _gff_type_iterator(feature, ['mRNA', 'transcript', 'inferred_parent']):
                logger.debug("Adding gene: {0}".format(gene))
                exons = []
                logger.debug("subfeatures: {0}".format(gene.sub_features))
                for exon in [f for f in gene.sub_features if f.type == 'exon']:
                    #link[gene.id] = link.setdefault(gene.id, 0) + 1
                    start = int(exon.location.start.position)# - 1    
                    end = int(exon.location.end.position)
                    strand = smap[exon.strand]
                    exons.append([chrom, start, end, strand])
                if len(exons) >= min_exons:
                    transcripts.append([gene.id, fname, exons])

    return transcripts

def read_bed_transcripts(bedfile, min_exons=1):
    names = {}
    transcripts = []
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
                    
                logger.debug("read_bed: adding {0}".format(vals[3]))
                exons = [[vals[0], 
                          chromStart + start, 
                          chromStart + start + size, 
                         vals[5]]
                         for start, size in zip(starts, sizes)
                         ]
                
                if len(exons) >= min_exons:
                    transcripts.append([name, bedfile, exons])
    return transcripts



def read_gff(fname, format="gtf", collection=None, prefix=None):
    #limits = dict(gff_type = ["mRNA", "exon"])
    smap = {"1":"+",1:"+","-1":"-",-1:"-"}
    link = {}
    for rec in GFF.parse(open(fname)):
        chrom = rec.id
        for feature in rec.features:
            for gene in _gff_type_iterator(feature, ['mRNA', 'transcript']):
                #sys.stderr.write("Adding gene: {0}\n".format(gene))
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
                    
                #sys.stderr.write("read_bed: adding {0}\n".format(vals[3]))
                for start, size in zip(starts, sizes):
                    #print     
                    #print vals[0], chromStart + start, chromStart + start + size, vals[5], name
                    if size < 0:
                        print line
                    logger.debug("read_bed: exon {0}:{1}-{2}".format(vals[0], chromStart + start, chromStart + start + size))
                    collection.add(vals[0], chromStart + start, chromStart + start + size, vals[5], "bla", transcript=name)    

    #print link
    for name in names.keys():
        collection.link_transcript_exons(name)
    return collection


from BCBio import GFF
import sys
import os
import logging 
from tempfile import NamedTemporaryFile
import subprocess as sp
import pysam
import pybedtools

def _create_tabix(fname, ftype):
    logger = logging.getLogger("pita")
    tabix_file = ""
    logger.info("Creating tabix index for %s", os.path.basename(fname))
    logger.debug("Preparing %s for tabix", fname)
    tmp = NamedTemporaryFile(prefix="pita", delete=False)
    preset = "gff"
    if ftype == "bed":
        cmd = "sort -k1,1 -k2g,2 {0} | grep -v track | grep -v \"^#\" > {1}"
        preset = "bed"
    elif ftype in ["gff", "gff3", "gtf"]:
        cmd = "sort -k1,1 -k4g,4 {0} | grep -v \"^#\" > {1}"

    # Sort the input file
    logger.debug(cmd.format(fname, tmp.name))
    sp.call(cmd.format(fname, tmp.name), shell=True)
    # Compress using bgzip
    logger.debug("compressing %s", tmp.name)
    tabix_file = tmp.name + ".gz"
    pysam.tabix_compress(tmp.name, tabix_file)
    tmp.close()
    # Index (using tabix command line, as pysam.index results in a Segmentation fault
    logger.debug("indexing %s", tabix_file)
    sp.call("tabix {0} -p {1}".format(tabix_file, preset), shell=True)
    return tabix_file

def exons_to_tabix_bed(exons):
    logger = logging.getLogger("pita")
    logger.debug("Converting %s exons to tabix bed", len(exons))
    tmp = NamedTemporaryFile(prefix="pita", delete=False)
    logger.debug("Temp name %s", tmp.name)
    for exon in sorted(exons, cmp=lambda x,y: cmp([x.chrom, x.start], [y.chrom, y.start])):
        tmp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            exon.chrom, exon.start, exon.end, exon.id, "0", exon.strand))

    tmp.close()
    tabix_fname = _create_tabix(tmp.name, "bed")
    return tabix_fname

def tabix_overlap(fname1, fname2, chrom, fraction):
    logger = logging.getLogger("pita")
    logger.debug("TABIX overlap between %s and %s, %s", fname1, fname2, fraction)

    tab1 = pysam.Tabixfile(fname1)
    tab2 = pysam.Tabixfile(fname2)

    if not ((chrom in tab1.contigs) and (chrom in tab2.contigs)):
        return
    
    fobj1 = TabixIteratorAsFile(tab1.fetch(chrom))
    tmp1 = NamedTemporaryFile(prefix="pita.", delete=False)
    for line in fobj1.readlines():
        tmp1.write("{}\n".format(line.strip()))
    
    fobj2 = TabixIteratorAsFile(tab2.fetch(chrom))
    tmp2 = NamedTemporaryFile(prefix="pita.", delete=False)
    for line in fobj2.readlines():
        tmp2.write("{}\n".format(line.strip()))

    tmp1.flush()
    tmp2.flush()

    b1 = pybedtools.BedTool(tmp1.name)
    b2 = pybedtools.BedTool(tmp2.name)

    intersect = b1.intersect(b2, f=fraction)

    tmp1.close()
    tmp2.close()
    
    for f in intersect:
        yield f

def merge_exons(starts, sizes, l=0):
    merge = []
    for i, (start1, start2, size) in enumerate(zip(starts[:-1], starts[1:], sizes[:-1])):
        if start1 + size + l >= start2:
            merge.append(i + 1)
    
    if len(merge) == 0:
        return starts, sizes
   
    i = len(starts) - 1
    new_starts = []
    new_sizes = []
    
    while i >= 0:
        if i in merge:
            j = i - 1
            while j in merge:
                j -= 1
            new_starts.append(starts[j])
            new_sizes.append(starts[i] - starts[j] + sizes[i])
            
            i = j - 1
        else:
            new_starts.append(starts[i])
            new_sizes.append(sizes[i])
            i -= 1
    
    return new_starts[::-1], new_sizes[::-1]

def _gff_type_iterator(feature, ftypes):
    if feature.type in ftypes:
        yield feature
    else:
        for feature in feature.sub_features:
            for f in _gff_type_iterator(feature, ftypes):
                yield f

def read_gff_transcripts(fobj, fname="", min_exons=1, merge=0):
    
    # Setup logging
    logger = logging.getLogger('pita')
  
    if merge > 0:
        logger.warning("Merging exons not yet implemented for GFF files!")

    #limits = dict(gff_type = ["mRNA", "exon"])
    smap = {"1":"+",1:"+","-1":"-",-1:"-", None:"+"}
    transcripts = []
    for rec in GFF.parse(fobj):
        chrom = rec.id
        for feature in rec.features:
            #logger.debug("feature: {0}", feature)
            
            for gene in _gff_type_iterator(feature, ['mRNA', 'transcript', 'inferred_parent']):
                #logger.debug("Adding gene: {0}", gene)
                exons = []
                #logger.debug("subfeatures: {0}", gene.sub_features)
                for exon in [f for f in gene.sub_features if f.type == 'exon']:
                    #link[gene.id] = link.setdefault(gene.id, 0) + 1
                    start = int(exon.location.start.position)# - 1    
                    end = int(exon.location.end.position)
                    strand = smap[exon.strand]
                    exons.append([chrom, start, end, strand])
                logger.debug("%s: %s - %s exons", fname, gene.id, len(exons))
                if len(exons) >= min_exons:
                    transcripts.append([gene.id, fname, exons])

    return transcripts

def read_bed_transcripts(fobj, fname="", min_exons=1, merge=0):
    
    # Setup logging
    logger = logging.getLogger('pita')
    
    names = {}
    transcripts = []
    line = fobj.readline()
    while line:
        if not line.startswith("track"):
            #logger.debug(line)
            try:
                vals = line.strip().split("\t")
                
                i = 1
                name = "%s|%s|%s" % (vals[0], vals[3], i)
                while names.has_key(name):
                    i += 1
                    name = "%s|%s|%s" % (vals[0], vals[3], i)
                names[name] = 1
                
                chromStart = int(vals[1])
                
                sizes = [int(x) for x in vals[10].strip(",").split(",")]
                starts = [int(x) for x in vals[11].strip(",").split(",")] 
                    
                starts, sizes = merge_exons(starts, sizes, l=merge)
                
                exons = [[vals[0], 
                          chromStart + start, 
                          chromStart + start + size, 
                          vals[5]]
                          for start, size in zip(starts, sizes)
                         ]
                    
                if len(exons) >= min_exons:
                    logger.debug("read_bed: adding %s", vals[3])
                    transcripts.append([name, fname, exons])
                else:
                    logger.debug("read_bed: not adding %s, filter on minimum exons", vals[3])
            
            except:
                print "Error parsing BED file"
                print line
                raise
                sys.exit(1)
        
        line = fobj.readline()
    return transcripts

class TabixIteratorAsFile:
    def __init__(self, x):
        self.x = x

    def read(self):
        return None

    def readline(self):
        try:
            return self.x.next()
        except StopIteration:
            return None

    def readlines(self):
        line = self.readline()
        while line:
            yield line
            line = self.readline()





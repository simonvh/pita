from pita.collection import Collection
from pita.io import TabixIteratorAsFile, read_gff_transcripts, read_bed_transcripts
import logging
import pysam
import itertools
import sys
from tempfile import NamedTemporaryFile

def get_chrom_models(chrom, anno_files, data, weight):
    sep = ":::"
    
    logger = logging.getLogger("pita")
    
    # Read annotation files
    mc = Collection()
    for name, fname, ftype in anno_files:
        logger.debug("Reading annotation from {0}".format(fname))
        tabixfile = pysam.Tabixfile(fname)
        #tabixfile = fname
        if chrom in tabixfile.contigs:
            fobj = TabixIteratorAsFile(tabixfile.fetch(chrom))
            if ftype == "bed":
                it = read_bed_transcripts(fobj, fname, 2)
            elif ftype in ["gff", "gtf", "gff3"]:
                it = read_gff_transcripts(fobj, fname, 2)
            for tname, source, exons in it:
                mc.add_transcript("{0}{1}{2}".format(name, sep, tname), source, exons)
            del fobj    
        tabixfile.close()
        del tabixfile

    for name, fname, span, extend in data:
        if span == "splice":
            logger.debug("Reading splice data {0} from {1}".format(name, fname))
            mc.get_splice_statistics(fname, name=name)
        else:
            logger.debug("Reading BAM data {0} from {1}".format(name, fname))
            mc.get_read_statistics(fname, name=name, span=span, extend=extend)
    
    models = {}
    exons = []
    for cluster in mc.get_connected_models():
        if len(cluster) == 0:
            continue
        best_model = mc.max_weight(cluster, weight)
        genename = "{0}:{1}-{2}_".format(
                                        best_model[0].chrom,
                                        best_model[0].start,
                                        best_model[-1].end,
                                        )
        
        rpkm = mc.get_weight(best_model, "RNAseq", "rpkm") 
        if rpkm >= 1:
            genename += "V"
        else:
            genename += "X"
       
        rpkm = mc.get_weight(best_model, "H3K4me3", "first")
        logger.debug("{0}: H3K4me3: {1}".format(genename, rpkm)) 
        if rpkm >= 1:
            genename += "V"
        else:
            genename += "X"

        other_exons = [e for e in set(itertools.chain.from_iterable(cluster)) if not e in best_model]  
        del cluster
         
        ### Ugly logging stuff
        best_ev = {}
        other_ev = {}
        for e in best_model:
            for ev in set([x.split(sep)[0] for x in e.evidence]):
                best_ev[ev] = best_ev.setdefault(ev, 0) + 1

        # Fast way to collapse
        for e in other_exons:
            for ev in set([x.split(sep)[0] for x in e.evidence]):
                other_ev[ev] = other_ev.setdefault(ev, 0) + 1
        ev = []
        for e in best_model + other_exons:
            for evidence in e.evidence:
                ev.append(evidence.split(sep))

        ### End ugly logging stuff
        logger.debug("Best model: {0} with {1} exons".format(genename, len(best_model)))
        models[genename] = [genename, best_model, ev, best_ev, other_ev]
        
        for exon in best_model:
            exons[exon] = genename
   
    
    del mc
    
    discard = {}
    overlap = get_overlapping_models(exons.keys())
    for e1, e2 in overlap:
        if not(discard.has_key(exons[e1]) or discard.has_key(exons[e2])):
            m1 = models[exons[e1]][1]
            m2 = models[exons[e2]][1]
            w1 = mc.get_weight(m1, "RNAseq", "total_rpkm")
            w2 = mc.get_weight(m2, "RNAseq", "total_rpkm")
            if w1 >= w2:
                discard[exons[e1]] = 1
            else:
                discard[exons[e2]] = 1
    
    print discard


    return models 
    

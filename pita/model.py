from pita.dbcollection import DbCollection
from pita.annotationdb import AnnotationDb
from pita.io import TabixIteratorAsFile, read_gff_transcripts, read_bed_transcripts
from pita.util import get_overlapping_models,to_genomic_orf,longest_orf
import logging
import pysam
import itertools
import sys
from tempfile import NamedTemporaryFile
from pita.config import SEP

def get_chrom_models(chrom, anno_files, data, weight, prune=None, index=None):
    
    logger = logging.getLogger("pita")
    
    try:
        # Read annotation files
        db = AnnotationDb(new=False, index=index)
        logger.info("Reading annotation for {0}".format(chrom))
        for name, fname, ftype, min_exons in anno_files:
            logger.info("Reading annotation from {0}".format(fname))
            tabixfile = pysam.Tabixfile(fname)
            #tabixfile = fname
            if chrom in tabixfile.contigs:
                fobj = TabixIteratorAsFile(tabixfile.fetch(chrom))
                if ftype == "bed":
                    it = read_bed_transcripts(fobj, fname, min_exons=min_exons, merge=10)
                elif ftype in ["gff", "gtf", "gff3"]:
                    it = read_gff_transcripts(fobj, fname, min_exons=min_exons, merge=10)
                for tname, source, exons in it:
                    db.add_transcript("{0}{1}{2}".format(name, SEP, tname), source, exons)
                del fobj    
            tabixfile.close()
            del tabixfile
        
        # Prune spurious exon linkages
        #for p in mc.prune():
        #    logger.debug("Pruning {0}:{1}-{2}".format(*p))

        mc = DbCollection(db, chrom)
        # Remove long exons with only one evidence source
        mc.filter_long(l=2000)
        # Remove short introns
        #mc.filter_short_introns()
        logger.info("Loading data for {0}".format(chrom))

        for name, fname, span, extend in data:
            if span == "splice":
                logger.info("Reading splice data {0} from {1}".format(name, fname))
                db.get_splice_statistics(fname, name=name)
            else:
                logger.info("Reading BAM data {0} from {1}".format(name, fname))
                db.get_read_statistics(fname, name=name, span=span, extend=extend, nreads=None)
        
        models = {}
        exons = {}
        logger.info("Calling transcripts for {0}".format(chrom))
        for cluster in mc.get_connected_models():
            while len(cluster) > 0:
                #logger.debug("best model")
                best_model = mc.max_weight(cluster, weight)
                #logger.debug("best variant")
                best_model = mc.get_best_variant(best_model, weight)                
                genename = "{0}:{1}-{2}_".format(
                                            best_model[0].chrom,
                                            best_model[0].start,
                                            best_model[-1].end,
                                            )
                   
                logger.info("Best model: {0} with {1} exons".format(genename, len(best_model)))
                models[genename] = [genename, best_model]
#            
#                for exon in best_model:
#                    exons[str(exon)] = [exon, genename]
#       
                for i in range(len(cluster) - 1, -1, -1):
                    if cluster[i][0].start <= best_model[-1].end and cluster[i][-1].end >= best_model[0].start:
                        del cluster[i]    
        discard = {}
#        if prune:
#            #logger.debug("Prune: {0}".format(prune))
#            overlap = get_overlapping_models([x[0] for x in exons.values()])
#            #logger.debug("{0} overlapping exons".format(len(overlap)))
#            
#            gene_count = {}
#            for e1, e2 in overlap:
#                gene1 = exons[str(e1)][1]
#                gene2 = exons[str(e2)][1]
#                gene_count[gene1] = gene_count.setdefault(gene1, 0) + 1
#                gene_count[gene2] = gene_count.setdefault(gene2, 0) + 1
#
#            for e1, e2 in overlap:
#                gene1 = exons[str(e1)][1]
#                gene2 = exons[str(e2)][1]
#                if not(discard.has_key(gene1) or discard.has_key(gene2)):
#                    m1 = models[gene1][1]
#                    m2 = models[gene2][1]
#                
#                    loc1,loc2 = sorted([m1, m2], cmp=lambda x,y: cmp(x[0].start, y[0].start))
#                    l1 = float(loc1[-1].end - loc1[0].start)
#                    l2 = float(loc2[-1].end - loc2[0].start)
#                    if loc2[-1].end > loc1[-1].end:
#                        overlap = float(loc1[-1].end - loc2[0].start)
#                    else:
#                        overlap = l2
#
#                    #logger.info("Pruning {} vs. {}".format(str(m1),str(m2)))
#                    logger.info("1: {}, 2: {}, overlap: {}".format(
#                        l1, l2, overlap))
#                    logger.info("Gene {} count {}, gene {} count {}".format(
#                        str(gene1), gene_count[gene1], str(gene2), gene_count[gene2]
#                        ))
#                   
#                    prune_overlap = 0.1
#                    if overlap / l1 < prune_overlap and overlap / l2 < prune_overlap:
#                        logger.info("Not pruning!")
#                        continue
#                    
#                    w1 = 0.0
#                    w2 = 0.0
#                    for d in prune:
#                        logger.debug("Pruning overlap: {0}".format(d))
#                        tmp_w1 = mc.get_weight(m1, d["name"], d["type"])
#                        tmp_w2 = mc.get_weight(m2, d["name"], d["type"])
#                        m = max((tmp_w1, tmp_w2))
#                        if m > 0:
#                            w1 += tmp_w1 / max((tmp_w1, tmp_w2))
#                            w2 += tmp_w2 / max((tmp_w1, tmp_w2))
#
#                    if w1 >= w2:
#                        discard[gene2] = 1
#                    else:
#                        discard[gene1] = 1
#        
#        del c
#    
        
        logger.info("Done calling transcripts for {0}".format(chrom))
        return [v for m,v in models.items() if not m in discard]

    except:
        logger.exception("Error on {0}".format(chrom))
  
    return []

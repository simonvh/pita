from pita.dbcollection import DbCollection
from pita.annotationdb import AnnotationDb
from pita.io import TabixIteratorAsFile, read_gff_transcripts, read_bed_transcripts
from pita.util import get_overlapping_models
import logging
import pysam
from pita.config import SEP

def load_chrom_data(conn, new, chrom, anno_files, data, index=None):
    logger = logging.getLogger("pita")
    
    try:
        # Read annotation files
        db = AnnotationDb(index=index, conn=conn, new=new)
        logger.debug("%s %s", chrom, id(db.session))
        logger.info("Reading annotation for %s", chrom)
        for name, fname, tabix_file, ftype, min_exons in anno_files:
            logger.info("Reading annotation from %s", fname)
            tabixfile = pysam.Tabixfile(tabix_file)
            #tabixfile = fname
            if chrom in tabixfile.contigs:
                fobj = TabixIteratorAsFile(tabixfile.fetch(chrom))
                if ftype == "bed":
                    it = read_bed_transcripts(fobj, fname, 
                            min_exons=min_exons, merge=10)
                elif ftype in ["gff", "gtf", "gff3"]:
                    it = read_gff_transcripts(fobj, fname, 
                            min_exons=min_exons, merge=10)
                for tname, source, exons in it:
                    db.add_transcript(
                            "{0}{1}{2}".format(name, SEP, tname), 
                            source, exons)
                del fobj    
            tabixfile.close()
            del tabixfile
        
        logger.info("Loading data for %s", chrom)

        for name, fname, span, extend in data:
            if span == "splice":
                logger.info("Reading splice data %s from %s", name, fname)
                db.get_splice_statistics(chrom, fname, name)
            else:
                logger.info("Reading BAM data %s from %s", name, fname)
                db.get_read_statistics(chrom, fname, name=name, span=span, extend=extend, nreads=None)
 
    except:
        logger.exception("Error on %s", chrom)
        raise

def get_chrom_models(conn, chrom, weight, repeats=None, prune=None, keep=None, filter_ev=None, experimental=None):
    if keep is None:
        keep = []
    if filter_ev is None:
        filter_ev = []
    if experimental is None:
        experimental = []

    logger = logging.getLogger("pita")
    logger.debug(str(weight)) 
    try:
        db = AnnotationDb(conn=conn)
        
        # Filter repeats
        if repeats:
            for x in repeats:
                db.filter_repeats(chrom, x)

        for ev in filter_ev:
            db.filter_evidence(chrom, ev, experimental) 
        
        mc = DbCollection(db, weight, prune=prune, chrom=chrom)
       
        # Remove short introns
        #mc.filter_short_introns()
      
        models = {}
        exons = {}
        logger.info("Calling transcripts for %s", chrom)
        for model in mc.get_best_variants(weight):
            genename = "{0}:{1}-{2}_".format(
                                        model[0].chrom,
                                        model[0].start,
                                        model[-1].end,
                                        )
                   
                
            logger.info("Best model: %s with %s exons", 
                    genename, len(model))
            models[genename] = [genename, model]
            
            for exon in model:
                exons[str(exon)] = [exon, genename]

        discard = {}
        if prune:
            logger.debug("Prune: {0}".format(prune))
            overlap = get_overlapping_models([x[0] for x in exons.values()])
            if len(overlap) > 1:
                logger.info("%s overlapping exons", len(overlap))
#                logger.warn("Overlap: {0}".format(overlap))
                
            gene_count = {}
            for e1, e2 in overlap:
                gene1 = exons[str(e1)][1]
                gene2 = exons[str(e2)][1]
                gene_count[gene1] = gene_count.setdefault(gene1, 0) + 1
                gene_count[gene2] = gene_count.setdefault(gene2, 0) + 1

            for e1, e2 in overlap:
                gene1 = exons[str(e1)][1]
                gene2 = exons[str(e2)][1]
                if not(gene1 in discard or gene2 in discard):
                    m1 = models[gene1][1]
                    m2 = models[gene2][1]
                
                    loc1,loc2 = sorted([m1, m2], cmp=lambda x,y: cmp(x[0].start, y[0].start))
                    l1 = float(loc1[-1].end - loc1[0].start)
                    l2 = float(loc2[-1].end - loc2[0].start)
                    if loc2[-1].end > loc1[-1].end:
                        overlap = float(loc1[-1].end - loc2[0].start)
                    else:
                        overlap = l2

                    #logger.info("Pruning {} vs. {}".format(str(m1),str(m2)))
                    #logger.info("1: {}, 2: {}, overlap: {}".format(
                    #    l1, l2, overlap))
                    #logger.info("Gene {} count {}, gene {} count {}".format(
                    #    str(gene1), gene_count[gene1], str(gene2), gene_count[gene2]
                    #    ))
#                   
                    prune_overlap = prune["overlap"]["fraction"]
                    if overlap / l1 < prune_overlap and overlap / l2 < prune_overlap:
                        logger.debug("Not pruning because fraction of overlap is too small!")
                        continue
                    
                    w1 = 0.0
                    w2 = 0.0
                    for d in prune["overlap"]["weights"]:
                        logger.debug("Pruning overlap: %s", d)
                        tmp_w1 = -mc.get_weight(m1)
                        tmp_w2 = -mc.get_weight(m2)
                        m = max((tmp_w1, tmp_w2))
                        if m > 0:
                            w1 += tmp_w1 / max((tmp_w1, tmp_w2))
                            w2 += tmp_w2 / max((tmp_w1, tmp_w2))

                    if w1 >= w2:
                        logger.info("Discarding %s", gene2)
                        discard[gene2] = 1
                    else:
                        logger.info("Discarding %s", gene1)
                        discard[gene1] = 1
        
        logger.info("Done calling transcripts for %s", chrom)
        result = [v for m,v in models.items() if not m in discard]
        #print "VV", result
        return [[name, [e.to_flat_exon() for e in exons]] for name, exons in result]

    except:
        logger.exception("Error on %s", chrom)
  
    return []

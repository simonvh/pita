from pita.collection import Collection
from pita.io import TabixIteratorAsFile, read_gff_transcripts, read_bed_transcripts
from pita.util import get_overlapping_models,to_genomic_orf,longest_orf
import logging
import pysam
import itertools
import sys
from tempfile import NamedTemporaryFile

def get_chrom_models(chrom, anno_files, data, weight, prune=None, index=None):
    sep = ":::"
    
    logger = logging.getLogger("pita")
    
    try:
        # Read annotation files
        mc = Collection(index)
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
        
        # Prune spurious exon linkages
        #for p in mc.prune():
        #    logger.debug("Pruning {0}:{1}-{2}".format(*p))

        # Remove long exons with only one evidence source
        mc.filter_long(l=2000)
        # Remove short introns
        mc.filter_short_introns()

        for name, fname, span, extend in data:
            if span == "splice":
                logger.debug("Reading splice data {0} from {1}".format(name, fname))
                mc.get_splice_statistics(fname, name=name)
            else:
                logger.debug("Reading BAM data {0} from {1}".format(name, fname))
                mc.get_read_statistics(fname, name=name, span=span, extend=extend, nreads=None)
        
        models = {}
        exons = {}
        for cluster in mc.get_connected_models():
            while len(cluster) > 0:
                best_model = mc.max_weight(cluster, weight)
                variants = [m for m in mc.all_simple_paths(best_model[0], best_model[-1])]
                if len(variants) > 1:
                    logger.info("Checking {0} extra variants".format(len(variants)))
                    best_model = mc.max_weight(variants, weight)
                
                 
                genename = "{0}:{1}-{2}_".format(
                                            best_model[0].chrom,
                                            best_model[0].start,
                                            best_model[-1].end,
                                            )
            
                rpkm = mc.get_weight(best_model, "RNAseq", "rpkm") 
                if rpkm >= 0.2:
                    genename += "V"
                else:
                    genename += "X"
           
                rpkm = mc.get_weight(best_model, "H3K4me3", "first_rpkm")
                logger.debug("{0}: H3K4me3: {1}".format(genename, rpkm)) 
                if rpkm >= 1:
                    genename += "V"
                else:
                    genename += "X"
    
                other_exons = [e for e in set(itertools.chain.from_iterable(cluster)) if not e in best_model]  
                for i in range(len(cluster) - 1, -1, -1):
                    if cluster[i][0].start <= best_model[-1].end and cluster[i][-1].end >= best_model[0].start:
                        del cluster[i]
                    
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
                    exons[str(exon)] = [exon, genename]
       
        
        discard = {}
        if prune:
            logger.debug("Prune: {0}".format(prune))
            overlap = get_overlapping_models([x[0] for x in exons.values()])
            logger.debug("{0} overlapping exons".format(len(overlap)))
            for e1, e2 in overlap:
                gene1 = exons[str(e1)][1]
                gene2 = exons[str(e2)][1]
                if not(discard.has_key(gene1) or discard.has_key(gene2)):
                    m1 = models[gene1][1]
                    m2 = models[gene2][1]
                
                    w1 = 0.0
                    w2 = 0.0
                    for d in prune:
                        logger.debug("Pruning overlap: {0}".format(d))
                        tmp_w1 = mc.get_weight(m1, d["name"], d["type"])
                        tmp_w2 = mc.get_weight(m2, d["name"], d["type"])
                        m = max((tmp_w1, tmp_w2))
                        if m > 0:
                            w1 += tmp_w1 / max((tmp_w1, tmp_w2))
                            w2 += tmp_w2 / max((tmp_w1, tmp_w2))

                    if w1 >= w2:
                        discard[gene2] = 1
                    else:
                        discard[gene1] = 1
        
        del mc
    
        return [v for m,v in models.items() if not m in discard]

    except:
        logger.exception("Error on {0}".format(chrom))
  
    return []

from pita.collection import Collection
from pita.io import TabixIteratorAsFile, read_gff_transcripts, read_bed_transcripts
import logging
import pysam


def get_chrom_models(chrom, anno_files, data, weight):
    logger = logging.getLogger("pita")
    
    # Read annotation files
    mc = Collection()
    for name, fname, ftype in anno_files:
        logger.debug("Reading annotation from {0}".format(fname))
        tabixfile = pysam.Tabixfile(fname)
        if chrom in tabixfile.contigs:
            fobj = TabixIteratorAsFile(tabixfile.fetch(chrom))
            if ftype == "bed":
                it = read_bed_transcripts(fobj, fname, 3)
            elif ftype in ["gff", "gtf", "gff3"]:
                it = read_gff_transcripts(fobj, fname, 3)
            for tname, source, exons in it:
                mc.add_transcript("{0}:{1}".format(name, tname), source, exons)

    for name, fname, span, extend in data:
        logger.debug("Reading data {0} from {1}".format(name, fname))
        mc.get_read_statistics(fname, name=name, span=span, extend=extend)

    for i, cluster in enumerate(mc.get_all_transcript_clusters()):
        best_model = mc.max_weight(cluster, weight)
        genename = "{0}:{1}-{2}_".format(
                                        best_model[0].chr,
                                        best_model[0].start,
                                        best_model[-1].end,
                                        )
        
        if mc.get_weight(best_model, "H3K4me3", "first") > 100:
            genename += "V"
        else:
            genename += "X"

        abs_x = mc.get_weight(best_model, "RNAseq", "all") 
        x = abs_x / (sum([e.end - e.start for e in best_model]) / 1000.0)
        if x > 5:
            genename += "V"
        else:
            genename += "X"
          
        yield genename, best_model, cluster


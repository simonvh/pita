from itertools import izip
import logging
import random
import re

import numpy as np
import networkx as nx 
from networkx.algorithms.connectivity import minimum_st_node_cut
from networkx.algorithms.flow import edmonds_karp

from pita.util import longest_orf,exons_to_seq

class DbCollection(object):
    def __init__(self, db, chrom=None):
        # dict with chrom as key
        self.logger = logging.getLogger("pita")
        
        self.db = db
        self.chrom = chrom

        # All transcript models will be stored as a directed (acyclic) graph
        self.graph = nx.DiGraph()

        # Store read counts of BAM files
        self.nreads = {}
        
        # Store extension used in BAM statistics
        self.extend = {}

        self.logger.debug("Loading exons in graph")
        for exon in self.db.get_exons(chrom):
            self.add_feature(exon)
        
        self.logger.debug("Loading introns in graph")
        n = 0
        #for junction in self.db.get_splice_junctions(chrom, ev_count=1, read_count=20):
        for junction in self.db.get_splice_junctions(chrom, ev_count=0, read_count=1):
            n += 1
            self.add_feature(junction)
        self.logger.debug("%s introns were loaded", n)

    def add_feature(self, feature):
        """ 
        """

        if feature.ftype == "exon":
            # Add chromosome to keys
            nodes = [
                (feature.in_node(), "exon_in"), 
                (feature.out_node(), "exon_out"),
                ]
            for node,ftype in nodes:
                f = self.graph.add_node(node, ftype=ftype)
            self.graph.add_path([n[0] for n in nodes], 
                    weight=-1, ftype='exon')

        elif feature.ftype == "splice_junction":
          
            # Add transcript model to the graph
            for e1,e2 in self.db.get_junction_exons(feature):
                if e1.strand == "+":
                    self.graph.add_path((e1.out_node(), e2.in_node()), 
                            ftype="splice_junction", weight=-1)
                if e2.strand == "-":
                    self.graph.add_path((e2.out_node(), e1.in_node()), 
                            ftype="splice_junction", weight=-1)
  
    def get_best_variants(self, weights):
        self._set_weights(weights)
        
        add = {}
        add_ends = []
        for i,model in enumerate(nx.weakly_connected_components(self.graph)):
            source = "source_{}".format(i + 1)
            sink = "sink_{}".format(i + 1)
            ends = [k for k,v in self.graph.out_degree(model).items() if v == 0]
            for end in ends:
                add_ends.append((end, sink))
            
            for node in model:
                if self.graph.node[node]['ftype'] ==  "exon_in":
                    add.setdefault(source, []).append(node)
        for p in add_ends:
            self.graph.add_path(p)
        
        for source, targets in add.items():
            for target in targets:
                self.graph.add_path((source, target), ftype='source')

            sink = source.replace("source", "sink")
            pred,dis = nx.bellman_ford(self.graph, source)
            t = sink
            best_variant = []
            while pred[t]:
                best_variant.append(pred[t])
                t = pred[t]
            
            p = re.compile(r'(.+):(\d+)([+-])')
            model = []
            for i in range(0, len(best_variant) - 1, 2):
                n1,n2 = best_variant[i:i+2]
                m = p.search(n1)
                chrom, start, strand = [m.group(x) for x in (1,2,3)]
                m = p.search(n2)
                start = int(start)
                end = int(m.group(2))
                if start > end:
                    end,start = start,end
                #chrom,start,_,strand = n1.split
                e = self.db.fetch_feature((chrom, start, end, strand, "exon", None))
                model.append(e)
            if strand == "+":
                model = model[::-1]
            yield model

    def _set_weights(self, weights):
      
        p = re.compile(r'(.+):(\d+)([+-])')
        for n1,n2 in self.graph.edges():
            d = self.graph.edge[n1][n2]
            if d['ftype'] == "exon":
                m = p.search(n1)
                chrom, start, strand = [m.group(x) for x in (1,2,3)]
                m = p.search(n2)
                start = int(start)
                end = int(m.group(2))
                if start > end:
                    end,start = start,end
                #chrom,start,_,strand = n1.split
                e = self.db.fetch_feature((chrom, start, end, strand, "exon", None))
                
                for iw in weights:
                    weight = iw["weight"]
                    idtype = iw["type"]
                    identifier = iw["name"]
    
                    signal = self.db.feature_stats(e, identifier)
                    length = end - start
                    if idtype == "all":
                        d[identifier] = signal
                    elif idtype == "rpkm":
                        if signal == 0:
                            d[identifier] = 0
                        else:
                            print signal, self.db.nreads(identifier)
                            d[identifier] = float(signal) / (self.db.nreads(identifier) / 1e6) / length * 1000.0 

            elif d['ftype'] == "source":
                


                
                print "SOURUCE"
                for iw in weights:
                    weight = iw["weight"]
                    idtype = iw["type"]
                    identifier = iw["name"]
                    if idtype == "first":
                        m = p.search(n2)
                        strand = m.group(3)
                        print "f", n2
                        #if strand == "+":
                        #    return signal[0]
                        #else:
                        #return signal[-1]
 

        print "flop start"
        for n1,n2 in self.graph.edges():
            d = self.graph.edge[n1][n2]
            if 'rnaSeq' in d:
                d['weight'] = -d['rnaSeq']
            print d
        print "flop end"
        return 
              

    def get_weight(self, transcript, identifier, idtype):
        signal = [self.db.feature_stats(e, identifier) for e in transcript]
        exon_lengths = [e.end - e.start for e in transcript]

        if idtype in ["all", "rpkm", "weighted"]:
            total_signal = float(sum(signal))
            total_exon_length = sum(exon_lengths)
            
            if idtype == "all":
                return total_signal
            elif idtype == "rpkm":
                if total_signal == 0:
                    return 0
                return float(total_signal) / (self.db.nreads(identifier) / 1e6) / total_exon_length * 1000.0 

            elif idtype == "weighted":
                return float(total_signal) / total_exon_length * len(transcript)

        elif idtype in ["mean_exon", "total_rpkm"]:
            all_exons = [s/float(l) for s,l in zip(signal, exon_lengths)]
            nreads = self.db.nreads(identifier)
            if not nreads:
                nreads = 1000000
                self.logger.warn("Number of reads in db is 0 for %s", identifier)
            rpkms = [s * 1000.0 / nreads * 1e6 for s in all_exons]
            if idtype == "mean_exon":
                if len(rpkms) == 0:
                    self.logger.warning("Empty score array for mean_exon")
                    return 0
                else:
                    return np.mean(rpkms)
 
            if idtype == "total_rpkm":
                return sum(rpkms) 
        
        elif idtype == "first":
            if transcript[0].strand == "+":
                return signal[0]
            else:
                return signal[-1]
       
        elif idtype == "first_rpkm":
            exon = transcript[0]
            count = signal[0]
            if transcript[0].strand == "-":
                exon = transcript[-1]
                count = signal[-1]

            size = exon.end - exon.start
            extend = self.extend.setdefault(identifier, (0,0))
            size += extend[0] +  extend[1]
            if count == 0:
                return 0
            
            rpkm = count / (self.db.nreads(identifier)/ 1e6) / size * 1000.0
            
            return rpkm  
        
        elif idtype == "splice":
            if len(transcript) == 1:
                return 0
            w = []
            for e1, e2 in zip(transcript[:-1], transcript[1:]):
                w.append(self.db.splice_stats(e1, e2, identifier))
            if len(w) == 0:
                self.logger.warning("Empty score array for splice")
                return 0
            else:
                return np.sum(w)
        
        elif idtype == "orf":
            start, end = longest_orf(exons_to_seq(transcript))
            return end - start

        elif idtype == "evidence":
            #return 1
            evidences = [len(exon.evidences) for exon in transcript]
            if len(evidences) == 0:
                self.logger.warning("Empty score array for evidence")
                return 0
            else:
                return np.mean(evidences)
    
        elif idtype == "length":
            return np.sum([e.end - e.start for e in transcript])
        
        else:
            raise Exception, "Unknown idtype"

    def max_weight(self, transcripts, identifier_weight):
        max_transcripts = 10000
        if len(transcripts) > max_transcripts:
            self.logger.warn("More than %s transcripts, random sampling to a managable number", max_transcripts)
            transcripts = random.sample(transcripts, max_transcripts)
        
        if not identifier_weight or len(identifier_weight) == 0:
            w = [len(t) for t in transcripts]    
        else:
            w = np.array([0] * len(transcripts))
            pseudo = 1e-10
            for iw in identifier_weight:
                weight = iw["weight"]
                idtype = iw["type"]
                identifier = iw["name"]
                
                idw = []
                for transcript in transcripts:
                    tw = self.get_weight(transcript, identifier, idtype)
                    idw.append(pseudo + tw)
    
                idw = np.array(idw)
                idw = idw / max(idw) * weight
                w = w + idw
        self.db.clear_stats_cache() 
        return transcripts[np.argmax(w)]

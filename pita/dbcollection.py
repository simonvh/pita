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
        for exon in self.db.get_exons(chrom, eager=True):
            self.add_feature(exon)
        
        self.logger.debug("Loading introns in graph")
        n = 0
        #for junction in self.db.get_splice_junctions(chrom, ev_count=1, read_count=20):
        for junction in self.db.get_splice_junctions(chrom, ev_count=0, read_count=1, eager=True):
            n += 1
            self.add_feature(junction)
        self.logger.debug("%s introns were loaded", n)

        self.max_id_value = {}

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
            self._set_edge_weight(feature, nodes[0][0], nodes[1][0], {})
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
        
        iweight = {}
        for iw in weights:
            identifier = iw["name"]
            weight = iw["weight"]
            self.max_id_value[identifier] = 0
            print "HOEHA", identifier, weight
            iweight[identifier] = weight
        
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
                self._set_source_weight((source, target), weights)
                print "ha", self.graph.edge[source][target]

        
        for n1,n2 in self.graph.edges():
            self.graph.edge[n1][n2]['weight'] = -0.01
            d = self.graph.edge[n1][n2]
            
            print d.keys()
            print d.values()
            print iweight.keys()
            print iweight.values()
            
            for k,v in self.max_id_value.items():
                if v > 0:
                    print "Mwaha", k, v
                    
                    
                    if k in d:
                        w = d[k] / v * iweight[k]
                        self.graph.edge[n1][n2]['weight'] -= w 
        
        print "hola"
        for n1,n2 in self.graph.edges():
            d = self.graph.edge[n1][n2]
            print n1, n2, d['weight']
        print "mola"
        
        for source, targets in add.items():
            sink = source.replace("source", "sink")
            pred,dis = nx.bellman_ford(self.graph, source)
            t = sink
            best_variant = []
            while pred[t]:
                best_variant.append(pred[t])
                t = pred[t]
            
            p = re.compile(r'(.+):(\d+)([+-])')
            model = []
            strand = "+"
            for i in range(0, len(best_variant) - 1, 2):
                n1,n2 = best_variant[i:i+2]
                print "Getting {} {}".format(n1, n2)
                e = self._nodes_to_exon(n1, n2)
                strand = e.strand
                model.append(e)
            if strand == "+":
                model = model[::-1]
            print model
            yield model

    def _nodes_to_feature(self, n1, n2, feature): 
        p = re.compile(r'(.+):(\d+)([+-])')
        m = p.search(n1)
        chrom, start, strand = [m.group(x) for x in (1,2,3)]
        m = p.search(n2)
        start = int(start)
        end = int(m.group(2))
        if start > end:
            end,start = start,end
        e = self.db.fetch_feature((chrom, start, end, strand, feature, None))
        return e
    
    def _nodes_to_exon(self, n1, n2): 
        return self._nodes_to_feature(n1, n2, "exon")

    def _nodes_to_splice_junction(self, n1, n2): 
        return self._nodes_to_feature(n1, n2, "splice_junction")
    
    def _set_edge_weight(self, feature, n1, n2, weights):
        d = self.graph.edge[n1][n2]
        
        feature_stats = {}
        for f in feature.read_counts: 
            name = f.read_source.name
            feature_stats[name] = feature_stats.get(name, 0) + f.count
       
        for iw in weights:
            weight = iw["weight"]
            idtype = iw["type"]
            identifier = iw["name"]
            id_value = 0
            signal = feature_stats.get(identifier, 0)
            if d['ftype'] == "exon":
                length = feature.end - feature.start
                if idtype == "all":
                    id_value = signal
                elif idtype == "rpkm":
                    if signal > 0:
                        mreads = self.db.nreads(identifier) / 1e6
                        id_value = float(signal) / mreads / length * 1000.0
                elif idtype == "evidence":
                    id_value = len (feature.evidences)
                elif idtype == "length":
                    id_value = feature.end - feature.start
            elif d['ftype'] == "splice_junction":
                if idtype == "splice":
                    f = self._nodes_to_splice_junction(n1, n2)
                    id_value = self.db.intron_splice_stats(f, identifier)
                    print "SPLICE:", id_value
                elif idtype == "orf":
                    raise NotImplementedError
                    #start, end = longest_orf(exons_to_seq(transcript))
                    #return end - start
            d[identifier] = id_value
            if id_value > self.max_id_value[identifier]:
                self.max_id_value[identifier] = id_value
        
    def _set_weights(self, weights):
      
        for n1,n2 in self.graph.edges():
            d = self.graph.edge[n1][n2]
            for iw in weights:
                weight = iw["weight"]
                idtype = iw["type"]
                identifier = iw["name"]
                id_value = 0
                if d['ftype'] == "exon":
                    e = self._nodes_to_exon(n1, n2) 
                    signal = self.db.feature_stats(e, identifier)
                    length = e.end - e.start
                    if idtype == "all":
                        id_value = signal
                    elif idtype == "rpkm":
                        if signal > 0:
                            mreads = self.db.nreads(identifier) / 1e6
                            id_value = float(signal) / mreads / length * 1000.0
                    elif idtype == "evidence":
                        id_value = len (e.evidences)
                    elif idtype == "length":
                        id_value = e.end - e.start
                elif d['ftype'] == "splice_junction":
                    if idtype == "splice":
                        f = self._nodes_to_splice_junction(n1, n2)
                        id_value = self.db.intron_splice_stats(f, identifier)
                        print "SPLICE:", id_value
                    elif idtype == "orf":
                        raise NotImplementedError
                        #start, end = longest_orf(exons_to_seq(transcript))
                        #return end - start
                d[identifier] = id_value
                if id_value > self.max_id_value[identifier]:
                    self.max_id_value[identifier] = id_value

    def _set_source_weight(self, edge, weights):
        n2 = edge[-1] 
        print "SOURCE"
        for iw in weights:
            weight = iw["weight"]
            idtype = iw["type"]
            identifier = iw["name"]
            if idtype == "first":
                print "f", n2
                for other in self.graph[n2]:
                    print "other", other 
                    e = self._nodes_to_exon(n2, other)
                    signal = self.db.feature_stats(e, identifier)
                    self.graph.edge[edge[0]][n2][identifier] = signal
                    if signal > self.max_id_value[identifier]:
                        self.max_id_value[identifier] = signal


    def get_weight(self, m):
        return 1
        #for e1, e2 in zip(m[:-1], m[1:]):
    
    def filter_long(self, l=1000, evidence=2):
        for exon in self.db.get_long_exons(self.chrom, l, evidence):
            out_edges = len(self.graph.out_edges([exon]))
            in_edges = len(self.graph.in_edges([exon]))
            self.logger.debug("Filter long: %s, in %s out %s", 
                    exon, in_edges, out_edges)

            if (in_edges >= 0 and out_edges >= 1 and exon.strand == "+" 
                    or in_edges >= 1 and out_edges >= 0 and exon.strand == "-"):
                self.logger.info("Removing long exon %s", exon)
                self.graph.remove_node(exon)


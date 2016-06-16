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
    def __init__(self, db, weights, prune=None, chrom=None):
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

        # store maximum weight per type
        self.max_id_value = {}
        for iw in weights:
            identifier = iw["name"]
            self.max_id_value[identifier] = 0
        
        # Load the exons
        self._load_exons(weights, chrom=chrom, prune=prune)
        
        # Load the introns
        self._load_splice_junctions(weights, chrom=chrom, prune=prune)

    def _load_exons(self, weights, chrom=None, prune=None):
        """
        Load exons from database in graph
        """
        self.logger.debug("Loading exons in graph")
        
        # Get filtering rules 
        if prune == None:
            prune = {}
        min_length = prune.get("exons", {}).get("min_length", None)
        max_length = prune.get("exons", {}).get("max_length", None)
        ev = prune.get("exons", {}).get("evidence", None)
        if min_length:
            self.logger.debug("Minimum length for exon: %s", 
                    min_length)
        if max_length:
            self.logger.debug("Maximum length for exon: %s", 
                    max_length)
        if ev:
            self.logger.debug("Minimum evidence for exon: %s", 
                    ev)
        
        for exon in self.db.get_exons(chrom, eager=True, 
                min_length=min_length, max_length=max_length, evidence=ev):
            self.add_feature(exon, weights)
    
    def _load_splice_junctions(self, weights, chrom=None, prune=None): 
        """
        Load splice junctions from database in graph
        """
        self.logger.debug("Loading introns in graph")
        
        # Get filtering rules 
        if prune == None:
            prune = {}
        min_reads = prune.get("introns", {}).get("min_reads", None)
        ev = prune.get("introns", {}).get("evidence", None)
        if min_reads:
            self.logger.debug("Minimum reads for splice junction: %s", 
                    min_reads)
        if ev:
            self.logger.debug("Minimum evidence for splice junction: %s", 
                    ev)
        
        # Load splice junctions    
        n = 0
        for junction in self.db.get_splice_junctions(chrom, 
                ev_count=ev, read_count=min_reads, eager=True):
            n += 1
            self.add_feature(junction, weights)
        self.logger.debug("%s introns were loaded", n)

    def add_feature(self, feature, weights):
        """ 
        Add feature to the graph
        """

        # Exon
        if feature.ftype == "exon":
            nodes = [
                    (feature.in_node(), "exon_in"), 
                    (feature.out_node(), "exon_out"),
                    ]
            for node,ftype in nodes:
                f = self.graph.add_node(node, ftype=ftype)
            self.graph.add_path([n[0] for n in nodes], 
                    weight=-1, ftype='exon')
            self._set_edge_weight(feature, nodes[0][0], nodes[1][0], weights)
        # Intron
        elif feature.ftype == "splice_junction":
            for e1,e2 in self.db.get_junction_exons(feature):
                if e1.strand == "-":
                    e1,e2 = e2,e1
                self.graph.add_path((e1.out_node(), e2.in_node()), 
                        ftype="splice_junction", weight=-1)
                self._set_edge_weight(feature, e1.out_node(), e2.in_node(), weights)

    def get_best_variants(self, weights):

        iweight = {}
        for iw in weights:
            identifier = iw["name"]
            weight = iw["weight"]
            iweight[identifier] = weight

        add = {}
        add_ends = []
        for i,model in enumerate(nx.weakly_connected_components(self.graph)):
            source = "source_{}".format(i + 1)
            sink = "sink_{}".format(i + 1)
            ends = [k for k,v in self.graph.out_degree(model).items() 
                    if self.graph.node[k].get('ftype', "") == "exon_out" and v == 0]
            for end in ends:
                add_ends.append((end, sink))

            for node in model:
                if self.graph.node[node].get('ftype', "") ==  "exon_in":
                    add.setdefault(source, []).append(node)

        for p in add_ends:
            self.graph.add_path(p, ftype="sink")

        for source, targets in add.items():
            for target in targets:
                self.graph.add_path((source, target), ftype='source')
                self._set_source_weight((source, target), weights)

        for n1,n2 in self.graph.edges():
            self.graph.edge[n1][n2]['weight'] = -0.01
            d = self.graph.edge[n1][n2]
            for k,v in self.max_id_value.items():
                if v > 0:
                    if k in d:
                        w = d[k] / float(v) * iweight[k]
                        self.graph.edge[n1][n2]['weight'] -= w 

            if self.graph.edge[n1][n2]['ftype'] == "exon":
                self.logger.debug("Edge: %s, %s", n1, n2)
                for k,v in d.items():
                    self.logger.debug("Key: %s, value: %s", k, v)

        for n1,n2 in self.graph.edges():
            d = self.graph.edge[n1][n2]

        for source, targets in add.items():
            sink = source.replace("source", "sink")
            try:
                pred,dis = nx.bellman_ford(self.graph, source)
                if not pred.has_key(sink):
                    continue
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
                    e = self._nodes_to_exon(n1, n2)
                    if e:
                        strand = e.strand
                        model.append(e)
                if strand == "+":
                    model = model[::-1]
                #print model
                if len(model) > 0:
                    yield model
    
            except Exception as e:
                raise
                self.logger.warning("Failed: %s", self.graph.edge[source].keys())
                self.logger.warning("%s", e)
    
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
                    id_value = length
                elif idtype == "orf":
                    start, end = longest_orf(feature.seq)
                    id_value =  end - start
            elif d['ftype'] == "splice_junction":
                if idtype == "splice":
                    f = self._nodes_to_splice_junction(n1, n2)
                    id_value = feature_stats.get(identifier, 0)
            d[identifier] = id_value
            if id_value > self.max_id_value[identifier]:
                self.max_id_value[identifier] = id_value
        
    def _set_source_weight(self, edge, weights):
        n2 = edge[-1] 
        for iw in weights:
            weight = iw["weight"]
            idtype = iw["type"]
            identifier = iw["name"]
            if idtype == "first":
                for other in self.graph[n2]:
                    e = self._nodes_to_exon(n2, other)
                    signal = self.db.feature_stats(e, identifier)
                    self.graph.edge[edge[0]][n2][identifier] = signal
                    if signal > self.max_id_value[identifier]:
                        self.max_id_value[identifier] = signal

    def _model_to_path(self, model):
        """
        takes a list of Exons and returns a list of edges
        """
        
        # path is reversed if on minus strand
        if model[0].strand == "-":
            model = model[::-1]
        
        # add source node to the path
        source = [n for n in self.graph.predecessors(model[0].in_node())
                    if n.startswith("source")][0]
        path = [
                self.graph.edge[source][model[0].in_node()],
                self.graph.edge[model[0].in_node()][model[0].out_node()],
                    ]
        # add edge for every exon and intron
        for e1,e2 in zip(model[:-1], model[1:]):
            path += [
                    self.graph.edge[e1.out_node()][e2.in_node()],
                    self.graph.edge[e2.in_node()][e2.out_node()]
                    ]
        return path

    def get_weight(self, m):
        """
        return the total weight for a model
        """

        return sum([edge.get('weight', 0) for edge in self._model_to_path(m)])
    
    def filter_long(self, l=1000, evidence=2):
        """
        remove exon with length of at least <l> that is not supported by at least
        <evidence> sources
        """
        
        for exon in self.db.get_long_exons(self.chrom, l, evidence):
            out_edges = len(self.graph.out_edges([exon]))
            in_edges = len(self.graph.in_edges([exon]))
            self.logger.debug("Filter long: %s, in %s out %s", 
                    exon, in_edges, out_edges)

            if (in_edges >= 0 and out_edges >= 1 and exon.strand == "+" 
                    or in_edges >= 1 and out_edges >= 0 and exon.strand == "-"):
                self.logger.info("Removing long exon %s", exon)
                self.graph.remove_node(exon)


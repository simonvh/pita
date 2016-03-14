from itertools import izip
import logging
import random

import numpy as np
import networkx as nx 
from networkx.algorithms.connectivity import minimum_st_node_cut
from networkx.algorithms.flow import edmonds_karp

from pita.util import longest_orf,exons_to_seq

def connected_models(graph):
    for u, v in graph.edges():
        graph[u][v]['weight'] = -1
    for c in nx.weakly_connected_components(graph):
        starts =  [k for k,v in graph.in_degree(c).items() if v == 0]
        ends = [k for k,v in graph.out_degree(c).items() if v == 0]
        paths = []
        
        for s in starts:
            order,d = nx.bellman_ford(graph,s, weight='weight')
            
            for e in ends:
                if e in d: 
                    path = [e]
                    x = e
                    while order[x]:
                        path.append(order[x])
                        x = order[x]
            
                    paths.append(path[::-1])
        yield paths

def recursive_neighbors(graph, node_list):
    result = []
    for node in node_list:
        result += graph.neighbors(node)
    new_list = list(set(node_list + result))
    for node in result:
        if node not in node_list:
            return recursive_neighbors(graph, new_list)
    return new_list


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
            self.graph.add_node(feature)
            self.graph.node[feature]['weight'] = 1
        
        elif feature.ftype == "splice_junction":
          
            # Add transcript model to the graph
            for exons in self.db.get_junction_exons(feature):
                self.graph.add_path(exons)
            
#    def remove_exon(self, e):
#        if e in self.graph:
#            self.logger.info("Removing exon {0}".format(e))
#            self.graph.remove_node(e)
#            del self.exons[e.chrom][to_sloc(e.start, e.end, e.strand)]

    def get_initial_exons(self, chrom=None):
        """ Return all leftmost exons
        """
        in_degree = self.graph.in_degree(self.get_exons(chrom)).items()
        return [k for k,v in in_degree if v == 0]

    def get_connected_models(self):
        for paths in connected_models(self.graph):
            if len(paths) > 0:
                self.logger.debug("yielding %s paths", len(paths))
            yield paths
   
    def get_node_cuts(self, model):
        node_cuts = []
        cuts = list(minimum_st_node_cut(self.graph, model[0], model[-1], flow_func=edmonds_karp))
        while len(cuts) == 1:
            node_cuts = cuts + node_cuts
            cuts = list(minimum_st_node_cut(self.graph, model[0], cuts[0], flow_func=edmonds_karp))
        return node_cuts

    def get_best_variant(self, model, weight):
        
        if len(model) == 1:
            return model

        nodeset = self.get_node_cuts(model)
        if len(list(nodeset)) > 0:
            self.logger.debug("option 1")
            self.logger.debug(str(nodeset))
            nodeset = [model[0]] + list(nodeset) + [model[-1]]
            self.logger.debug("got nodeset")
            best_variant = [model[0]]
            for n1,n2 in zip(nodeset[:-1], nodeset[1:]):
                self.logger.debug("%s %s", str(n1), str(n2))
                variants = [m for m in self.all_simple_paths(n1, n2)]
                self.logger.debug("Got %s variants", len(variants))
                best_variant += self.max_weight(variants, weight)[1:]
                self.logger.debug("Best variant: %s", best_variant)
        else:
            variants = [m for m in self.all_simple_paths(model[0], model[-1])]
            best_variant = self.max_weight(variants, weight)
      
        e = best_variant[-1]
        if e.strand == "+":
            best_variant[-1] = self.db.get_longest_3prime_exon(e.chrom, e.start, e.strand)
        else:
            e = best_variant[0]
            best_variant[0] = self.db.get_longest_3prime_exon(e.chrom, e.end, e.strand)
        return best_variant

    def prune(self):
        pruned = []

        for cluster in self.get_connected_models():
            self.logger.debug("Pruning %s models", len(cluster))
            
            discard = []
            new_cluster = [m for m in cluster]
            
            while len(new_cluster) > 0:
                #print len(new_cluster)
                #c_min = min([m[0].start for m in new_cluster])
                #c_max = max([m[-1].end for m in new_cluster])
                #print c_min, c_max
                #selection = [m for m in new_cluster if m[0].start == c_min or m[-1].end == c_max]
                
                longest = sorted(new_cluster, cmp=lambda x,y: cmp(x[-1].end - x[0].start, y[-1].end - y[0].start))[-1]
                discard.append(longest)
                new_cluster = [m for m in new_cluster if m != longest]
                if len(new_cluster) != 0:
                    
                    graph = nx.DiGraph()
                    for m in new_cluster:
                        graph.add_path(m)
                    
                    result = [x for x in connected_models(graph) if len([y for y in x if len(y) > 2]) > 1]
                    if len(result) > 1:
                        break
            
            if len(new_cluster) != 0:
                #print len(new_cluster)
                discard_edges = []
                for m in discard:
                    for e1, e2 in zip(m[:-1], m[1:]):
                        discard_edges.append((e1, e2))
                
                keep_edges = []
                for m in new_cluster:
                    for e1, e2 in zip(m[:-1], m[1:]):
                        keep_edges.append((e1, e2))

                for x in set(discard_edges) - set(keep_edges):
                    self.graph.remove_edge(x[0], x[1])

                    pruned.append([x[0].chrom, x[0].end, x[1].start])
        
        return pruned

    def is_weak_splice(self, splice, evidence=1):
        exons = []
        my = []
        for e1,e2 in self.db.get_junction_exons(splice):
            if e1 in self.graph and e2 in self.graph and (e1, e2) in self.graph.edges():
                my.append((e1.end, e2.start))        
                for e in [e1, e2]:
                    if e not in exons:
                        exons.append(e) 
        
        if len(exons) == 0:
            return False

        splices = self.graph.edges(recursive_neighbors(self.graph, exons))
        if len(splices) == 1:
            return False
        
        counts = {}
        for s in splices:
            if (s[0].end, s[1].start) not in counts:
                counts[(s[0].end, s[1].start)] = self.db.get_splice_count(s[0], s[1])

        bla = [v for k,v in counts.items() if k not in my]
        if len(bla) == 0:
            return False
        self.logger.debug("%s %s %s", counts[my[0]], np.mean(bla),  np.std(bla))
        return counts[my[0]] < 0.1 * np.mean(bla)# - np.std(bla))

    def prune_splice_junctions(self, max_reads=10, evidence=2, keep=None):
        if keep is None:
            keep = []

        keep = set(keep)
        for splice in self.db.get_splice_junctions(self.chrom, max_reads=max_reads):
            self.logger.debug("Splice %s, evidence %s", splice, len(splice.evidences))
            ev_sources = [e.source for e in splice.evidences]
            if keep.intersection(ev_sources):
                self.logger.debug("Keeping this splice %s", splice)
                continue
            if len(splice.evidences) <= evidence:
                self.logger.debug("Checking splice %s", splice)
                if self.is_weak_splice(splice, evidence):
                    self.logger.debug("Removing splice %s", splice)
                    for e1,e2 in self.db.get_junction_exons(splice):
                        if (e1,e2) in self.graph.edges():
                            self.graph.remove_edge(e1, e2)
                            for node in e1, e2:
                                if len(self.graph.edges(node)) == 0:
                                    self.logger.debug("Removing lonely exon %s", node)
                                    self.graph.remove_node(node)
   
    def filter_long(self, l=1000, evidence=2):
        for exon in self.db.get_long_exons(self.chrom, l, evidence):
            out_edges = len(self.graph.out_edges([exon]))
            in_edges = len(self.graph.in_edges([exon]))
            self.logger.debug("Filter long: %s, in %s out %s", exon, in_edges, out_edges)

            if in_edges >= 0 and out_edges >= 1 and exon.strand == "+" or in_edges >= 1 and out_edges >= 0 and exon.strand == "-":
                self.logger.info("Removing long exon %s", exon)
                self.graph.remove_node(exon)
    
    def filter_and_merge(self, nodes, l):
        for e1, e2 in self.graph.edges_iter(nodes):
            if e2.start - e1.end <= l:
                new_exon = self.add_exon(e1.chrom, e1.start, e2.end, e1.strand)
                self.logger.info("Adding %s", new_exon)
                for e_in in [e[0] for e in self.graph.in_edges([e1])]:
                    self.graph.add_edge(e_in, new_exon)
                for e_out in [e[1] for e in self.graph.out_edges([e2])]:
                    self.graph.add_edge(new_exon, e2)

                for e in (e1, e2):
                    self.remove_exon(e)                         
                
                return new_exon 
        return None

    def filter_short_introns(self, l=10, mode='merge'):
        filter_nodes = []
        for intron in self.graph.edges_iter():
            e1,e2 = intron
            if e2.start - e1.end <= l:
                filter_nodes += [e1, e2]
            
        if mode == "merge":
            exon = self.filter_and_merge(filter_nodes, l)
            while exon:
                exon = self.filter_and_merge(filter_nodes + [exon], l)
        else:
            for e in filter_nodes:
                self.remove_exon(e)                         
                    
    def all_simple_paths(self, exon1, exon2):
        return nx.all_simple_paths(self.graph, exon1, exon2)
    
    def get_alt_splicing_exons(self):
        for exon in self.get_exons():
            out_exons = [e[1] for e in self.graph.out_edges([exon]) if len(self.graph.out_edges([e[1]])) > 0]
            if len(out_exons) > 1:
                
                out_exon = out_exons[0]
                
                self.logger.info("ALT SPLICING %s %s", exon, out_exon)

            #in_exons = [e[0] for e in self.graph.in_edges([exon])]
            #for in_exon in in_exons:
            #    my_in_exons = [e[0] for e in self.graph.in_edges([in_exon])]
            #    for my_in_exon in my_in_exons:
            #        if my_in_exon in in_exons:
            #            self.logger.info("{0} is alternative exon".format(in_exon)) 

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

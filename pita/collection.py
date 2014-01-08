from pita.exon import *
from pita.util import read_statistics
import numpy
import sys
import logging
import pickle
from networkx.algorithms.components.connected import connected_components
import networkx as nx 
from itertools import izip, count

def to_loc(chrom, start, end, strand):
    return "{0}:{1}{3}{2}".format(chrom, start, end, strand)

def to_sloc(start, end, strand):
    return "{0}{2}{1}".format(start, end, strand)

def connected_models(graph):
    for u, v in graph.edges():
        graph[u][v]['weight'] = -1
    
    for c in nx.weakly_connected_components(graph):
        starts =  [k for k,v in graph.in_degree(c).items() if v == 0]
        ends = [k for k,v in graph.out_degree(c).items() if v == 0]
        paths = []
        
        for i,s in enumerate(starts):
            order,d = nx.bellman_ford(graph,s, weight='weight')
            
            for e in ends:
                if d.has_key(e): 
                    path = [e]
                    x = e
                    while order[x]:
                        path.append(order[x])
                        x = order[x]
            
                    paths.append(path[::-1])
        yield paths

class Collection:
    def __init__(self):
        # dict with chrom as key
        self.exons = {}
        self.logger = logging.getLogger("pita")
        
        # All transcript models will be stored as a directed (acyclic) graph
        self.graph = nx.DiGraph()

        # Store read counts of BAM files
        self.nreads = {}
        
        # Store extension used in BAM statistics
        self.extend = {}

    def add_exon(self, chrom, start, end, strand):
        """ 
        Create Exon object if does not exist.
        Returns Exon object.
        """

        # Add chromosome to keys
        self.exons.setdefault(chrom, {})

        # Return existing exon if it's there
        if self.exons[chrom].has_key(to_sloc(start, end, strand)):
            return self.exons[chrom][to_sloc(start, end, strand)]
        
        # Otherwise create new exon and return it
        e = Exon(chrom, start, end, strand)
        self.exons[chrom][to_sloc(start, end, strand)] = e
        return e

    def get_exons(self, chrom=None):
        """ Return exon objects
        Optional argument is the chromosome, otherwise return all
        """
        if chrom:
            return self.exons[chrom].values()
        else:
            exons = []
            for chrom in self.exons.keys():
                exons += self.exons[chrom].values()
            return exons

    def add_transcript(self, name, source, exons):
        """
        Add a transcript to the collection
        """

        # Sanity checks
        for e1, e2 in zip(exons[:-1], exons[1:]):
            if e1[0] != e2[0]:
                sys.stderr.write("{0} - {1}\n".format(e1, e2))
                raise ValueError, "Different chromosomes!"
            if e2[1] <= e1[2]:
                sys.stderr.write("{0} - {1}\n".format(e1, e2))
                raise ValueError, "exons overlap, or in wrong order"
            if e1[3] != e2[3]:
                sys.stderr.write("{0} - {1}\n".format(e1, e2))
                raise ValueError, "strands don't match"
        
        # First add all exons
        exons = [self.add_exon(*exon) for exon in exons]
        
        # Add transcript model to the graph
        self.graph.add_path(exons)
        
        # Default node weight = 1
        for e in exons:
            self.graph.node[e]['weight'] = 1

        [exon.add_evidence(name) for exon in exons]

    def get_initial_exons(self, chrom=None):
        """ Return all leftmost exons
        """
        in_degree = self.graph.in_degree(self.get_exons(chrom)).items()
        return [k for k,v in in_degree if v == 0]

    def get_connected_models(self):
        for paths in connected_models(self.graph):
            if len(paths) > 0:
                self.logger.debug("yielding {0} paths".format(len(paths)))
            yield paths
    
    def prune(self):
        pruned = []

        for i,cluster in enumerate(self.get_connected_models()):
            self.logger.debug("Pruning {0} models".format(len(cluster)))
            #print i + 1
            
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

    def get_read_statistics(self, fnames, name, span="exon", extend=(0,0)):
        from fluff.fluffio import get_binned_stats
        from tempfile import NamedTemporaryFile

        self.extend[name] = extend
        tmp = NamedTemporaryFile()
        estore = {}
        self.logger.debug("Writing exons to file")
        for exon in self.get_exons():
            start = exon.start
            end = exon.end
            if exon.strand == "-":
                start -= extend[1]
                end += extend[0]
            else:
                start -= extend[0]
                end += extend[1]
            if start < 0:
                start = 0
            
            estore["%s:%s-%s" % (exon.chrom, start, end)] = exon
            tmp.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                exon.chrom,
                start,
                end,
                str(exon),
                0,
                exon.strand
            ))
        tmp.flush()

        if type("") == type(fnames):
            fnames = [fnames]
       
        if not self.nreads.has_key(name):
            self.nreads[name] = 0

        for fname in fnames:
            if fname.endswith("bam"):
                rmrepeats = True
                self.nreads[name] += read_statistics(fname) 
            else:
                rmrepeats = False
            
            self.logger.debug("Getting overlap from {0}".format(fname))
            result = get_binned_stats(tmp.name, fname, 1, rpkm=False, rmdup=True, rmrepeats=rmrepeats)
        
            self.logger.debug("Reading results, save to exon stats")
        
            for row in result:
                vals = row.strip().split("\t")
                e = "%s:%s-%s" % (vals[0], vals[1], vals[2])
                c = float(vals[3])
                estore[e].stats[name] = estore[e].stats.setdefault(name, 0) + c
                self.graph.node[estore[e]][name] = -estore[e].stats[name]
            
        tmp.close()

    def get_splice_statistics(self, fnames, name):
        if type("") == type(fnames):
            fnames = [fnames]
        
        nrsplice = {}
        for fname in fnames:
            for line in open(fname):
                vals = line.strip().split("\t")
                chrom = vals[0]
                start,end,count = [int(x) for x in vals[1:]]
                
                if not nrsplice.has_key(chrom):
                    nrsplice[chrom] = {}
            
                if not nrsplice[chrom].has_key(start):
                    nrsplice[chrom][start] = {}
        
                nrsplice[chrom][start][end] =  nrsplice[chrom][start].setdefault(end, 0) +  count
        
        #print nrsplice[scaffold_1][
        for exon in self.get_exons():
            splice_start = exon.end
            exon.stats[name] = {}
            if nrsplice.has_key(exon.chrom):
                for end, count in nrsplice[exon.chrom].setdefault(splice_start, {}).items():
                    exon.stats[name][end] = count
                #print exon, exon.stats[name]

    def get_weight(self, transcript, identifier, idtype):
        if idtype == "all":
            total_signal = sum([e.stats.setdefault(identifier, 0) for e in transcript])
            return float(total_signal)
        
        elif idtype == "rpkm":
            total_exon_length = sum([e.end - e.start for e in transcript])
            total_signal = sum([e.stats.setdefault(identifier, 0) for e in transcript])
            return float(total_signal) / (self.nreads[identifier] / 1e6) / total_exon_length * 1000.0 

        elif idtype == "weighted":
            total_exon_length = sum([e.end - e.start for e in transcript])
            total_signal = sum([e.stats.setdefault(identifier, 0) for e in transcript])
            return float(total_signal) / total_exon_length * len(transcript)

        elif idtype == "mean_exon":
            all_exons = [e.stats.setdefault(identifier, 0) / (e.end - e.start) for e in transcript]
            return mean(all_exons)

        elif idtype == "total_rpkm":
            all_exons = [e.stats.setdefault(identifier, 0) / (e.end - e.start) for e in transcript]
            return sum(all_exons) 
        
        elif idtype == "splice":
            #self.logger.debug("SPLICE!")
            w = 0.0
            for e1, e2 in zip(transcript[:-1], transcript[1:]):
                w += e1.stats[identifier].setdefault(e2.start, 0)
            #self.logger.debug("Weight: {0}".format(w))
            return w

        elif idtype == "first":
            if transcript[0].strand == "+":
                return transcript[0].stats.setdefault(identifier,0)
            else:
                return transcript[-1].stats.setdefault(identifier,0)

        elif idtype == "first_rpkm":
            if transcript[0].strand == "+":
                exon = transcript[0]
            else:
                exon = transcript[-1]

            size = exon.end - exon.start
            size += self.extend[identifier][0] +  self.extend[identifier][1]
            count = exon.stats.setdefault(identifier, 0)
            rpkm = count / (self.nreads[identifier] / 1e6) / size * 1000.0
            
            return rpkm  
        

    def max_weight(self, transcripts, identifier_weight):
        #identifier_weight = []
        
        if not identifier_weight or len(identifier_weight) == 0:
            w = [len(t) for t in transcripts]    
            #sys.stderr.write("weights: {0}".format(w))
        else:
            w = numpy.array([0] * len(transcripts))
            pseudo = 1e-10
            for iw in identifier_weight:
                weight = iw["weight"]
                idtype = iw["type"]
                identifier = iw["name"]
                
                idw = []
                for transcript in transcripts:
                   tw = self.get_weight(transcript, identifier, idtype)
                   self.logger.debug("{0} - {1} - {2}".format(
                                                             transcript[0],
                                                             identifier,
                                                             tw
                                                             ))
                   idw.append(pseudo + tw)
    
                idw = numpy.array(idw)
                idw = idw / max(idw) * weight
                #sys.stderr.write("Adding {0}\n".format(idw))
                w = w + idw
        
        return transcripts[numpy.argmax(w)]

def get_updated_exons(model, name):
    strand = model[0].strand
    n = 0
    u5 = 0
    u3 = 0
    for e in model:
        if name in [x.split(":")[0] for x in e.evidence]:
            break
        n += 1

    if strand == "-":
        u3 = n
    else:
        u5 = n 
     
    n = 0
    for e in model[::-1]:
        if name in [x.split(":")[0] for x in e.evidence]:
            break
        n += 1

    if strand == "-":
        u5 = n
    else:
        u3 = n
    
    return u5,u3

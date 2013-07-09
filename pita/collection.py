from pita.exon import *
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

class Collection:
    def __init__(self):
        # dict with chrom as key
        self.exons = {}
        self.logger = logging.getLogger("pita")
        
        # All transcript models will be stored as a directed (acyclic) graph
        self.graph = nx.DiGraph()


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
                sys.stderr.write("{0} - {1}\n" % (e1, e2))
                raise ValueError, "Different chromosomes!"
            if e2[1] <= e2[2]:
                sys.stderr.write("{0} - {1}\n" % (e1, e2))
                raise ValueError, "exons overlap, or in wrong order"
            if e1[3] != e2[3]:
                sys.stderr.write("{0} - {1}\n" % (e1, e2))
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
        #self.logger.debug("get_connected_models")
        for u, v in self.graph.edges():
            self.graph[u][v]['weight'] = -1
        
        for c in nx.weakly_connected_components(self.graph):
            self.logger.debug("calculating paths of {0} exons".format(len(c)))
            starts =  [k for k,v in self.graph.in_degree(c).items() if v == 0]
            #self.logger.debug("{0} starts".format(len(starts)))
            ends = [k for k,v in self.graph.out_degree(c).items() if v == 0]
            #self.logger.debug("{0} ends".format(len(ends)))
            paths = []
            
            for i,s in enumerate(starts):
                #self.logger.debug("{0} out of {1} starts".format(i+ 1, len(starts)))
                #self.logger.debug("Starting at {0} ".format(str(s)))

                order,d = nx.bellman_ford(self.graph,s, weight='weight')
                
                for e in ends:
                    if d.has_key(e): 
                        path = [e]
                        x = e
                        while order[x]:
                            path.append(order[x])
                            x = order[x]
                
                        paths.append(path[::-1])

            self.logger.debug("yielding {0} paths".format(len(paths)))
            yield paths

    def get_read_statistics(self, fname, name, span="exon", extend=(0,0)):
        from fluff.fluffio import get_binned_stats
        from tempfile import NamedTemporaryFile

        tmp = NamedTemporaryFile()
        estore = {}
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

        if fname.endswith("bam"):
            rmrepeats = True
        else:
            rmrepeats = False
        
        result = get_binned_stats(tmp.name, fname, 1, rpkm=False, rmdup=True, rmrepeats=rmrepeats)
        
        for row in result:
            vals = row.strip().split("\t")
            e = "%s:%s-%s" % (vals[0], vals[1], vals[2])
            c = float(vals[3])
            estore[e].stats[name] = c
            self.graph.node[estore[e]][name] = -c

    def get_weight(self, transcript, identifier, idtype):
        if idtype == "all":
            return sum([exon.stats.setdefault(identifier, 0) for exon in transcript])
        elif idtype == "first":
            if transcript[0].strand == "+":
                return transcript[0].stats.setdefault(identifier,0)
            else:
                return transcript[-1].stats.setdefault(identifier,0)


    def max_weight(self, transcripts, identifier_weight):
        identifier_weight = []
        
        if len(identifier_weight) == 0:
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
                    idw.append(pseudo + self.get_weight(transcript, identifier, idtype))
    
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

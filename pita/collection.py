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
        
        # Create new exon and return it
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
        #sys.stderr.write("Adding {0} with {1} exons\n".format(name, len(exons)))
        # First add all exons
        exons = [self.add_exon(*exon) for exon in exons]
        
        self.graph.add_path(exons)
        for e in exons:
            self.graph.node[e]['lala'] = -1

        [exon.add_evidence(name) for exon in exons]

        # Now add evidence and link them
        for e1, e2 in zip(exons[0:-1], exons[1:]):
            try:
                e1.add_link(e2, name)
            except:
                self.logger.warn("Not linking {0} and {1}".format(str(e1), str(e2)))
    
    def get_initial_exons(self):
        """ Return all leftmost exons
        """
        return [exon for exon in self.get_exons() if not exon.linked]

    def get_paths(self, exon):
        return [exon.linked_chains(ev) for ev in exon.evidence]
    
    def get_transcripts(self):
        t = []
        for exon in self.get_initial_exons():
            for path in self.get_paths(exon):
                t += path
        
        return t

    def extend(self, paths):
        """Recursively extend transcript to all possible rightmost exons
        """
        new_paths = []
        if len([1 for path in paths if path[-1].get_linked_exons()]) > 0:
            for path in paths:
                exons = path[-1].get_linked_exons()
                if exons:
                    bla = [path + [exon] for exon in exons]
                    new_paths += self.extend(bla)
                else:
                    new_paths += [path]
            return new_paths
        else:
            return paths

    def get_all_transcripts(self):
        t = [] 
        #sys.stderr.write("Initial exons:\n")
        for exon in self.get_initial_exons():
            #sys.stderr.write("{0}\n".format(str(exon)))
            
            path = [[exon]]
            t += self.extend(path)
        
        #sys.stderr.write("****\n")
        return t

    def has_overlapping_exon(self, t1, t2):
        """ Return True if any exon in t1 is identical to any exon in t2
        """
        x = t1 + t2
        if len(set(t1 + t2)) != len(x):
            return True
        return False

    def get_overlapping_index(self, cluster, t):
        """Returns index of list of transcript cluster where an exon in t
        overlaps with an exon in any transcript of a cluster
        """
        for i, c in enumerate(cluster):
            for t2 in c:
                if self.has_overlapping_exon(t, t2):
                #    sys.stderr.write("Yep!\n")
                    return i

                #else:
                #    sys.stderr.write("No\n{0}\n{1}!\n".format(t,t2))
    
    def get_overlap_index(self, transcripts):
        m = 0
        for i, t1 in izip(count(), transcripts[:-1]):
            m = max([m, t1[-1].end])
            t2 = transcripts[i + 1]
            if t1[0].chr != t2[0].chr or t2[0].start > m:
                return i + 1
        return len(transcripts)
            
    def to_graph(self, l):
        G = nx.Graph()
        for part in l:
            G.add_nodes_from(part)
            G.add_edges_from(izip(part[0:-1], part[1:]))
        return G
   
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

                #for e in ends:
                #    for path in nx.all_simple_paths(self.graph, s, e):
                #        self.logger.debug("Adding {0}".format(str(path)))
                #        paths.append(path)
            
            self.logger.debug("yielding {0} paths".format(len(paths)))
            yield paths

    def get_all_transcript_clusters(self):
        """ Returns a list of lists of transcripts
        All transcript sharing at least one exon are grouped together
        """
       
        transcripts = self.get_all_transcripts()
        if len(transcripts ) > 0:
            self.logger.debug("Sorting transcripts ({0})".format(transcripts[0][0].chr))
            transcripts = sorted(transcripts, 
                             lambda x,y: cmp(
                                            (x[0].chr, x[0].start),
                                            (y[0].chr, y[0].start)
                                            )
                            )
            self.logger.debug("Done")
        
            self.logger.debug("Get overlapping transcripts")
        
            #f = open("transcripts.pickle", "w")
            #pickle.dump(transcripts, f)
            #f.close()
            
            #clusters = []
            while 1: 
                self.logger.debug("Overlap of {0}:{1}-{2}, {3} models".format(
                                                                transcripts[0][0].chr,
                                                                transcripts[0][0].start,
                                                                transcripts[1][-1].end,
                                                                len(transcripts)))
 
                overlap_i = self.get_overlap_index(transcripts)
                self.logger.debug("exon index")
                e_index = {}
                for i, t in izip(count(), transcripts[:overlap_i]):
                    for e in t:
                        e_index.setdefault(e, []).append(i)
                
                self.logger.debug("Building graph")
                G = self.to_graph(e_index.values())
                result = connected_components(G)
                for c in result:
                    yield [transcripts[x] for x in c]
                

                del transcripts[:overlap_i]
                if len(transcripts) == 0:
                    break
            self.logger.debug("Done")
        
           
        #return clusters 

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
            
            estore["%s:%s-%s" % (exon.chr, start, end)] = exon
            tmp.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                exon.chr,
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

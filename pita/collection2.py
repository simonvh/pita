from pita.exon import *
from numpy import argmax,array
import sys

def to_loc(chrom, start, end, strand):
    return "{0}:{1}{3}{2}".format(chrom, start, end, strand)

def to_sloc(start, end, strand):
    return "{0}{2}{1}".format(start, end, strand)

class Collection:
    def __init__(self):
        # dict with chrom as key
        self.exons = {}
    
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
        [exon.add_evidence(name) for exon in exons]

        # Now add evidence and link them
        for e1, e2 in zip(exons[0:-1], exons[1:]):
            try:
                e1.add_link(e2, name)
            except:
                sys.stderr.write("Not linking {0} and {1}\n".format(str(e1), str(e2)))
    
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
    
    def get_overlapping_cluster(self, transcripts):
        if len(transcripts) == 0:
            return [],[]
        c = [transcripts.pop(0)]
        while transcripts and transcripts[0][0].chr == c[-1][0].chr and transcripts[0][0].start < max([x[-1].end for x in c]):
            c.append(transcripts.pop(0))
        return c, transcripts

    
    def get_all_transcript_clusters(self):
        """ Returns a list of lists of transcripts
        All transcript sharing at least one exon are grouped together
        """
        clusters = []
        transcripts = self.get_all_transcripts()
        sys.stderr.write("Sorting\n")
        transcripts = sorted(transcripts, 
                             lambda x,y: cmp(
                                            (x[0].chr, x[0].start),
                                            (y[0].chr, y[0].start)
                                            )
                            )
        sys.stderr.write("Done\n")
        
        while 1: 
            cluster, transcripts = self.get_overlapping_cluster(transcripts)
            #print len(cluster), len(transcripts)
            tmp = []
            for t in cluster:
                i = self.get_overlapping_index(tmp, t)
                if i >= 0:
                    #print "Adding {0} to {1}".format(t, i)
                    tmp[i].append(t)
                else:
                    #print "Appending {0}".format(t)
                    tmp.append([t])
            clusters += tmp
            
            if len(transcripts) == 0:
                break
        
        
           
        return clusters 

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

    def get_weight(self, transcript, identifier, idtype):
        
        if idtype == "all":
            return sum([exon.stats.setdefault(identifier, 0) for exon in transcript])
        elif idtype == "first":
            if transcript[0].strand == "+":
                return transcript[0].stats.setdefault(identifier,0)
            else:
                return transcript[-1].stats.setdefault(identifier,0)


    def max_weight(self, transcripts, identifier_weight, idtype):
        w = array([0] * len(transcripts))
        pseudo = 1e-10
        for identifier,weight in identifier_weight.items():
            idw = []
            for transcript in transcripts:
                idw.append(pseudo + self.get_weight(transcript, identifier, idtype[identifier]))

            idw = array(idw)
            idw = idw / max(idw) * weight
            #sys.stderr.write("Adding {0}\n".format(idw))
            w = w + idw
        return transcripts[argmax(w)]

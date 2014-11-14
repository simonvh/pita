import os
import sys
import logging
from gimmemotifs.genome_index import GenomeIndex
from sqlalchemy import and_
from sqlalchemy.orm import joinedload
from pita import db_session
from pita.db_backend import * 
from pita.util import read_statistics, get_splice_score
import yaml

class AnnotationDb():
    def __init__(self, session=None, conn='mysql://pita:@localhost/pita', new=False, index=None):
        self.logger = logging.getLogger("pita")
        if session:
            self.session = session
        else:
            if conn.startswith("sqlite"):
                self.Session = db_session(conn, new)
                self.session = self.Session()
                self.engine = db_session.engine
            else:
                self._init_session(conn, new) 
        
        if index:
            self.index = GenomeIndex(index)
        else:
            self.index = None
    
    #def __destroy__(self):
    #    self.session.close()

    def _init_session(self, conn, new=False):
        self.engine = create_engine(conn)
        self.engine.raw_connection().connection.text_factory = str
        if new:
            Base.metadata.drop_all(self.engine)
            Base.metadata.create_all(self.engine)
        Base.metadata.bind =self.engine
        Session = scoped_session(sessionmaker(bind=self.engine))
        self.session = Session()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.session.close()
    
    def dump_yaml(self):
        dump_dict = {}
        dump_dict['feature'] = [[f.id, f.chrom.encode('ascii','ignore'), f.start, f.end, f.strand.encode('ascii','ignore'), f.ftype.encode('ascii','ignore'), f.seq.encode('ascii','ignore')] for f in self.session.query(Feature)]
        
        dump_dict['read_source'] = [[r.id, r.name.encode('ascii','ignore'), r.source.encode('ascii','ignore'), r.nreads] for r in self.session.query(ReadSource)]

        dump_dict['read_count'] = [[r.read_source_id, r.feature_id, r.count,r.span.encode('ascii','ignore'),r.extend_up, r.extend_down] for r in self.session.query(FeatureReadCount)]
        
        dump_dict['evidence'] = [[r.id, r.name.encode('ascii','ignore'), r.source.encode('ascii','ignore')] for r in self.session.query(Evidence)]
        
        dump_dict['feature_evidence'] = [[r.feature_id, r.evidence_id] for r in self.session.query(FeatureEvidence)]

        return yaml.dump(dump_dict)
    
    
    def load_yaml(self, fname):
        data = yaml.load(open(fname))
        source_map = {}
        for old_id,name,fname,nreads in data['read_source']:
            r = get_or_create(self.session, ReadSource,
                    name=name, source=fname, nreads=nreads)
            self.session.commit()
            source_map[old_id] = r.id
    
    
        t = ["chrom","start","end","strand","ftype","seq"]
        result = self.engine.execute(
            Feature.__table__.insert(),
            [dict(zip(t, row[1:])) for row in data['feature']]
            )
    
        first = self.fetch_feature(data['feature'][0][1:])
        last = self.fetch_feature(data['feature'][-1][1:])
    
        f_map = dict(zip([x[0] for x in data['feature']], range(first.id, last.id + 1)))
        data['read_count'] = [
                [source_map[row[0]]] + [f_map[row[1]]] + row[2:] for row in data['read_count']
            ]
        t = ["read_source_id", "feature_id", "count", "span", "extend_up", "extend_down"]
    
        result = self.engine.execute(
            FeatureReadCount.__table__.insert(),
            [dict(zip(t, row)) for row in data['read_count']]
            )
    
        t = ["name","source"]
        result = self.engine.execute(
            Evidence.__table__.insert(),
            [dict(zip(t, row[1:])) for row in data['evidence']]
            )
    
        first = self.fetch_evidence(data['evidence'][0][1:])
        last = self.fetch_evidence(data['evidence'][-1][1:])
    
        ev_map = dict(zip([x[0] for x in data['evidence']], range(first.id, last.id + 1)))
    
        data['feature_evidence'] = [
                [f_map[row[0]], ev_map[row[1]]] for row in data['feature_evidence']
            ]
    
        t = ["feature_id", "evidence_id"]
        result = self.engine.execute(
            FeatureEvidence.__table__.insert(),
            [dict(zip(t, row)) for row in data['feature_evidence']]
            )
    
    def add_transcript(self, name, source, exons):
        """
        Add a transcript to the database
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
       
        chrom = exons[0][0]
        strand = exons[0][-1]

        seqs = []
        for exon in exons:
            seq = ""
            if self.index:
                seq = ""
                if exon[1] - 20 > 0:
                    seq = self.index.get_sequence(chrom, exon[1] - 20, exon[2] + 20, strand)
                seqs.append(seq)
            exon = get_or_create(self.session, Feature,
                             chrom = chrom,
                             start = exon[1],
                             end = exon[2],
                             strand = strand,
                             ftype = "exon",
                             seq = seq[20:-20]
                             ) 
            exon.evidences.append(Evidence(name=name, source=source))

        splice_donors = []
        splice_acceptors = []
        for i,(start,end) in enumerate([(e1[2], e2[1]) for e1, e2 in zip(exons[0:-1], exons[1:])]):
            sj = get_or_create(self.session, Feature,
                             chrom = chrom,
                             start = start,
                             end = end,
                             strand = strand,
                             ftype = "splice_junction"
                             ) 
            sj.evidences.append(Evidence(name=name, source=source))
            
            if strand == "+":
                if len(seqs) > (i + 1) and len(seqs[i]) > 46:
                    splice_donors.append(["{}_{}".format(name, i + 1), seqs[i][-23:-14]])
                if len(seqs) > (i + 2) and len(seqs[i + 1]) > 46:
                    f = ["{}_{}".format(name, i + 1), seqs[i + 1][:23]]
                    splice_acceptors.append(f)
            else:
                if len(seqs) > (i + 2) and len(seqs[i + 1]) > 46:
                    f = ["{}_{}".format(name, i + 1), seqs[i + 1][-23:-14]]
                    splice_donors.append(f)
                     
                if len(seqs) > (i + 1) and len(seqs[i]) > 46:
                    f = ["{}_{}".format(name, i + 1), seqs[i][:23]]
                    splice_acceptors.append(f)
        
        donor_score = get_splice_score(splice_donors, 5)
        acceptor_score = get_splice_score(splice_acceptors, 3)
        if donor_score + acceptor_score < 0:
            self.logger.error("Skipping {}, splicing not OK!".format(name))
            self.session.rollback()
        else:
            self.session.commit()
    
    def get_features(self, ftype=None, chrom=None):
        self.session.query(Feature).options(
                joinedload('read_counts')).all()
        
        query = self.session.query(Feature)
        if chrom:
            query = query.filter(Feature.chrom == chrom)
        if ftype:
            query = query.filter(Feature.ftype == ftype)
        features = [f for f in query]
        return features

    def get_exons(self, chrom=None):
        return self.get_features(ftype="exon", chrom=chrom)
    
    def get_splice_junctions(self, chrom=None):
        return self.get_features(ftype="splice_junction", chrom=chrom)

    def get_long_exons(self, l):
        query = self.session.query(Feature)
        query = query.filter(Feature.ftype == 'exon')
        query = query.filter(Feature.end - Feature.start >= l)
        return [e for e in query if len(e.evidences) == 1]

#    def get_evidence_count(self, exon):
#        query = self.session.query(Feature)
#        query = query.filter(Feature.ftype == 'exon')
#        query = query.filter(Feature.end - Feature.start >= l)
#        return [e for e in query if len(e.evidences) == 1]
    
    def get_read_statistics(self, chrom, fnames, name, span="all", extend=(0,0), nreads=None):
        from fluff.fluffio import get_binned_stats
        from tempfile import NamedTemporaryFile

        if span not in ["all", "start", "end"]:
            raise Exception("Incorrect span: {}".format(span))
        
        tmp = NamedTemporaryFile()
        estore = {}
        self.logger.debug("Writing exons to file")
        for exon in self.get_exons(chrom):
            start = exon.start
            end = exon.end
            if span == "start":
                if exon.strand == "+":
                    end = start
                elif exon.strand == "-":
                    start = end
            if span == "end":
                if exon.strand == "+":
                    start = end
                elif exon.strand == "-":
                    end = start
            
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

        for i, fname in enumerate(fnames):
            
            read_source = get_or_create(self.session, ReadSource, name=name, source=fname)
            self.session.commit() 
            if fname.endswith("bam") and (not nreads or not nreads[i]):
                rmrepeats = True
                self.logger.debug("Counting reads in {0}".format(fname))
                read_source.nreads = read_statistics(fname)
            else:
                rmrepeats = False

            self.logger.debug("Getting overlap from {0}".format(fname))
            result = get_binned_stats(tmp.name, fname, 1, rpkm=False, rmdup=False, rmrepeats=False)

            self.logger.debug("Reading results, save to exon stats")

            for row in result:
                vals = row.strip().split("\t")
                e = "%s:%s-%s" % (vals[0], vals[1], vals[2])
                c = float(vals[3])
                exon = estore[e]
                
                count = get_or_create(self.session, FeatureReadCount,
                            feature_id = exon.id,
                            read_source_id = read_source.id,
                            span = span,
                            extend_up = extend[0],
                            extend_down = extend[1])
                self.session.commit()
                if not count.count:
                    count.count = c
                else:
                    count.count += c

            self.session.commit()
        tmp.close()

    def get_splice_statistics(self, chrom, fnames, name):
        if type("") == type(fnames):
            fnames = [fnames]

        nrsplice = {}
        for fname in fnames:
            read_source = get_or_create(self.session, ReadSource, name=name, source=fname)
            self.session.commit()
            for line in open(fname):
                vals = line.strip().split("\t")
                if vals[0] == chrom:
                    start, end, c = [int(x) for x in vals[1:4]]
                    strand = vals[5]
                    
                    splice = get_or_create(self.session, Feature,
                             chrom = chrom,
                             start = start,
                             end = end,
                             strand = strand,
                             ftype = "splice_junction"
                             ) 
                    self.session.commit()

                    count = get_or_create(self.session, FeatureReadCount,
                            feature_id = splice.id,
                            read_source_id = read_source.id)
                
                    if not count.count:
                        count.count = c
                    else:
                        count.count += c
            
                    self.session.commit()    
    
    def get_junction_exons(self, junction):
        
        left = self.session.query(Feature).filter(and_(
            Feature.chrom == junction.chrom,
            Feature.strand == junction.strand,
            Feature.end == junction.start
            ))
        
        right = self.session.query(Feature).filter(and_(
            Feature.chrom == junction.chrom,
            Feature.strand == junction.strand,
            Feature.start == junction.end
            ))

        exon_pairs = []
        for e1 in left:
            for e2 in right:
                exon_pairs.append((e1, e2))
        return exon_pairs

    def feature_stats(self, feature, identifier):
        q = self.session.query(FeatureReadCount, ReadSource).join(ReadSource)
        q = q.filter(FeatureReadCount.feature_id == feature.id)
        q = q.filter(ReadSource.name == identifier)
       
        return sum([x[0].count for x in q.all()])

    def splice_stats(self, exon1, exon2, identifier):
        q = self.session.query(Feature)
        q = q.filter(Feature.ftype == "splice_junction")
        q = q.filter(Feature.chrom == exon1.chrom)
        q = q.filter(Feature.strand == exon1.strand)
        q = q.filter(Feature.start == exon1.end)
        q = q.filter(Feature.end == exon2.start)

        splice = q.first()
        return self.feature_stats(splice, identifier)

    def nreads(self, identifier):
        q = self.session.query(ReadSource)
        q = q.filter(ReadSource.name == identifier)
        return sum([s.nreads for s in q.all()])

    def fetch_feature(self, f):
        """ Feature as list """
        chrom, start, end, strand, ftype, seq = f
        feature = self.session.query(Feature).filter(and_(
                Feature.chrom == chrom,
                Feature.start == start,
                Feature.end == end,
                Feature.strand == strand,
                Feature.ftype == ftype,
                Feature.seq == seq,
                ))
    
        return feature[0]
    
    def fetch_evidence(self, f):
        """ Feature as list """
        name, source = f
        evidence = self.session.query(Evidence).filter(and_(
            Evidence.name == name,
            Evidence.source == source
            ))

        return evidence[0]

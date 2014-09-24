import os
import sys
import logging
from gimmemotifs.genome_index import GenomeIndex
from sqlalchemy import Column, ForeignKey, Integer, String, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship,sessionmaker
from sqlalchemy import create_engine, and_
from sqlalchemy.ext.associationproxy import association_proxy

Base = declarative_base()

class FeatureEvidence(Base):
    __tablename__ = 'feature_evidence'
    feature_id = Column(Integer, ForeignKey('feature.id'), primary_key=True)
    evidence_id = Column(Integer, ForeignKey('evidence.id'), primary_key=True)
    evidence = relationship("Evidence")
    feature = relationship("Feature")

class Feature(Base):
    __tablename__ = 'feature'
    __table_args__ = (
             UniqueConstraint(
                'chrom', 
                'start', 
                'end', 
                'strand', 
                'ftype',
                name='uix_1'),
            )
    id = Column(Integer, primary_key=True)
    chrom = Column(String(250), nullable=False) 
    start = Column(Integer, nullable=False) 
    end = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)
    ftype = Column(String(250), nullable=False) 
    _evidences = relationship('FeatureEvidence')
    evidences = association_proxy('_evidences', 'evidence',
                    creator=lambda _i: FeatureEvidence(evidence=_i),
                )

class Evidence(Base):
    __tablename__ = "evidence"
    id = Column(Integer, primary_key=True)
    name = Column(String(50))
    source = Column(String(50))

class ReadSource(Base):
    __tablename__ = "read_source"
    id = Column(Integer, primary_key=True)
    name = Column(String(250))
    source = Column(String(250))

class FeatureReadcount(Base):
    __tablename__ = "read_count"
    id = Column(Integer, primary_key=True)
    read_source_id = Column(Integer, ForeignKey('read_source.id'), primary_key=True)
    feature_id = Column(Integer, ForeignKey('feature.id'), primary_key=True)
    read_source = relationship("ReadSource")
    feature = relationship("Feature")
    count = Column(Integer)
    span = Column(String(50))
    extend_up = Column(Integer)
    extend_down = Column(Integer)

#class Read(Base):
#    __tablename__ = 'read'
#    id = Column(Integer, primary_key=True)
#    name = Column(String(250), nullable=False)

#class SplicedRead(Base):
#    __tablename__ = 'spliced_read'
#    id = Column(Integer, primary_key=True)
#    splice_id = Column(Integer, ForeignKey('splice.id'))
#    read_id = Column(Integer, ForeignKey('read.id'))

def get_or_create(session, model, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance
    else:
        instance = model(**kwargs)
        session.add(instance)
        return instance

class AnnotationDb():
    def __init__(self, new=False):
        self.logger = logging.getLogger("pita")

        self.engine = create_engine('sqlite:///pita_test.db')
        if new:
            Base.metadata.drop_all(self.engine)
        Base.metadata.create_all(self.engine)
        Base.metadata.bind = self.engine
        DBSession = sessionmaker(bind=self.engine)
        self.session = DBSession()

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

        for exon in exons:
            exon = get_or_create(self.session, Feature,
                             chrom = chrom,
                             start = exon[1],
                             end = exon[2],
                             strand = strand,
                             ftype = "exon"
                             ) 
            exon.evidences.append(Evidence(name=name, source=source))

        for start,end in [(e1[2], e2[1]) for e1, e2 in zip(exons[0:-1], exons[1:])]:
            sj = get_or_create(self.session, Feature,
                             chrom = chrom,
                             start = start,
                             end = end,
                             strand = strand,
                             ftype = "splice_junction"
                             ) 
            sj.evidences.append(Evidence(name=name, source=source))
        self.session.commit()

    def get_exons(self, chrom=None):
        if chrom:
            it = self.session.query(Feature).filter(and_(Feature.chrom == chrom, Feature.ftype == "exon"))
        else:
            it = self.session.query(Feature).filter(Feature.ftype == "exon")
        exons = [e for e in it]
        return exons

    def get_read_statistics(self, fnames, name, span="exon", extend=(0,0), nreads=None):
        from fluff.fluffio import get_binned_stats
        from tempfile import NamedTemporaryFile

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

        #if not self.nreads.has_key(name):
        #    self.nreads[name] = 0

        for i, fname in enumerate(fnames):
            if fname.endswith("bam") and (not nreads or not nreads[i]):
                rmrepeats = True
                self.logger.debug("Counting reads in {0}".format(fname))
                #self.nreads[name] += read_statistics(fname)
            else:
                rmrepeats = False

            self.logger.debug("Getting overlap from {0}".format(fname))
            result = get_binned_stats(tmp.name, fname, 1, rpkm=False, rmdup=False, rmrepeats=False)

            self.logger.debug("Reading results, save to exon stats")

            for row in result:
                vals = row.strip().split("\t")
                print vals
            #    e = "%s:%s-%s" % (vals[0], vals[1], vals[2])
            #    c = float(vals[3])
            #    estore[e].stats[name] = estore[e].stats.setdefault(name, 0) + c
            #    self.graph.node[estore[e]][name] = -estore[e].stats[name]

        tmp.close()


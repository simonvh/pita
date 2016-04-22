from warnings import warn
from sqlalchemy import Column, ForeignKey, Integer, String, Text, UniqueConstraint, Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship,sessionmaker,scoped_session
from sqlalchemy import and_, event
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method

from pita.exon import Exon


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
    seq = Column(Text(), default="") 
    flag = Column(Boolean(), default=False)
    _evidences = relationship('FeatureEvidence')
    evidences = association_proxy('_evidences', 'evidence',
                    creator=lambda _i: FeatureEvidence(evidence=_i),
                )
    read_counts = relationship('FeatureReadCount')
   
    @hybrid_property
    def length(self):
        return self.end - self.start

    def __str__(self):
        return self.to_loc()
    
    def to_loc(self):
        return "{}:{}{}{}".format(
                self.chrom, self.start, self.strand, self.end
                )
    
    def to_sloc(self):
        return "{}{}{}".format(
                self.start, self.strand, self.end
                )
    
    def in_node(self):
        if self.strand == "+":
            return "{}:{}{}_{}".format(
                    self.chrom, self.start, self.strand, self.ftype + "_in")
        else: 
            return "{}:{}{}_{}".format(
                    self.chrom, self.end, self.strand, self.ftype + "_in")
    
    def out_node(self):
        if self.strand == "-":
            return "{}:{}{}_{}".format(
                    self.chrom, self.start, self.strand, self.ftype + "_out")
        else: 
            return "{}:{}{}_{}".format(
                    self.chrom, self.end, self.strand, self.ftype + "_out")
  

        
    def to_flat_exon(self):
        e = Exon(self.chrom, self.start, self.end, self.strand)
        e.seq = self.seq
        return e

    def overlap(self, exon, strand=True, fraction=False):
        if fraction:
            warn("the 'fraction' argument in overlap is not yet implemented")

        if strand and self.strand != exon.strand:
            return 0
        if exon.start >= self.start and exon.start <= self.end:
            return self.end - exon.start
        if exon.end >= self.start and exon.end <= self.end:
            return exon.end - self.start
        return 0

class Evidence(Base):
    __tablename__ = "evidence"
    id = Column(Integer, primary_key=True)
    name = Column(String(250))
    source = Column(String(250))

class ReadSource(Base):
    __tablename__ = "read_source"
    id = Column(Integer, primary_key=True)
    name = Column(String(250))
    source = Column(String(250))
    nreads = Column(Integer, default=0)

class FeatureReadCount(Base):
    __tablename__ = "read_count"
    read_source_id = Column(Integer, ForeignKey('read_source.id'), primary_key=True)
    feature_id = Column(Integer, ForeignKey('feature.id'), primary_key=True)
    read_source = relationship("ReadSource")
    feature = relationship("Feature")
    count = Column(Integer, default=0)
    span = Column(String(50), default="all", primary_key=True)
    extend_up = Column(Integer, default=0, primary_key=True)
    extend_down = Column(Integer, default=0, primary_key=True)

def get_or_create(session, model, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance
    else:
        instance = model(**kwargs)
        session.add(instance)
        return instance

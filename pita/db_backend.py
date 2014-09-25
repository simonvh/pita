from sqlalchemy import Column, ForeignKey, Integer, String, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship,sessionmaker,mapper
from sqlalchemy import create_engine, and_, event
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.inspection import inspect

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
    read_counts = relationship('FeatureReadCount')
    
    #read_counts = association_proxy('_read_counts', 'read_count',
    #                creator=lambda _i: FeatureReadCount(read_count=_i),
    #            )

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

class FeatureReadCount(Base):
    __tablename__ = "read_count"
    read_source_id = Column(Integer, ForeignKey('read_source.id'), primary_key=True)
    feature_id = Column(Integer, ForeignKey('feature.id'), primary_key=True)
    read_source = relationship("ReadSource")
    feature = relationship("Feature")
    count = Column(Integer, default=0)
    span = Column(String(50), default="all")
    extend_up = Column(Integer, default=0)
    extend_down = Column(Integer, default=0)

def get_or_create(session, model, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance
    else:
        instance = model(**kwargs)
        session.add(instance)
        return instance

import os
import sys
import logging
from gimmemotifs.genome_index import GenomeIndex
from sqlalchemy import or_,and_,func
from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session,sessionmaker,subqueryload
from pita.db_backend import Base,get_or_create,ReadSource,Feature,\
        FeatureReadCount,Evidence,FeatureEvidence 
from pita.util import read_statistics, get_splice_score
import yaml
from pita.io import exons_to_tabix_bed, tabix_overlap
from fluff.fluffio import get_binned_stats
from tempfile import NamedTemporaryFile

class AnnotationDb(object):
    def __init__(self, session=None, conn='mysql://pita:@localhost/pita', new=False, index=None):
        self.logger = logging.getLogger("pita")
        
        # initialize db session
        if session:
            self.session = session
        else:
            self._init_session(conn, new) 
        
        # index to retrieve sequence
        self.index = None
        if index:
            self.index = GenomeIndex(index)
   
        self.cache_splice_stats = {}
        self.cache_feature_stats = {}

    def _init_session(self, conn, new=False):
        self.engine = create_engine(conn)
        self.engine.raw_connection().connection.text_factory = str
        
        # recreate database
        if new:
            Base.metadata.drop_all(self.engine)
            Base.metadata.create_all(self.engine)
        
        Base.metadata.bind =self.engine
        Session = scoped_session(sessionmaker(bind=self.engine))
        self.session = Session()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.session.close()
    
    def dump_yaml(self):
        dump_dict = {}
        dump_dict['feature'] = [
                    [f.id, 
                    f.chrom.encode('ascii','ignore'), 
                    f.start, 
                    f.end, 
                    f.strand.encode('ascii','ignore'), 
                    f.ftype.encode('ascii','ignore'), 
                    f.seq.encode('ascii','ignore')] 
                for f in self.session.query(Feature)]
        
        dump_dict['read_source'] = [
                [r.id, 
                    r.name.encode('ascii','ignore'), 
                    r.source.encode('ascii','ignore'), 
                    r.nreads] 
                for r in self.session.query(ReadSource)]

        dump_dict['read_count'] = [
                [r.read_source_id, 
                    r.feature_id, 
                    r.count,r.span.encode('ascii','ignore'),
                    r.extend_up, 
                    r.extend_down]
                for r in self.session.query(FeatureReadCount)]
        
        dump_dict['evidence'] = [
                [r.id, r.name.encode('ascii','ignore'), 
                    r.source.encode('ascii','ignore')] 
                for r in self.session.query(Evidence)]
        
        dump_dict['feature_evidence'] = [
                [r.feature_id, r.evidence_id] 
                for r in self.session.query(FeatureEvidence)]

        return yaml.dump(dump_dict)
    
    def load_yaml(self, fname):
        data = yaml.load(open(fname))
        source_map = {}
        
        if not data['feature']:
            return 

        for old_id,name,fname,nreads in data['read_source']:
            r = get_or_create(self.session, ReadSource,
                    name=name, source=fname, nreads=nreads)
            self.session.commit()
            source_map[old_id] = r.id
    
    
        t = ["chrom","start","end","strand","ftype","seq"]
        self.engine.execute(
            Feature.__table__.insert(),
            [dict(zip(t, row[1:])) for row in data['feature']]
            )
         
        self.session.commit()
        
        first = self.fetch_feature(data['feature'][0][1:])
        last = self.fetch_feature(data['feature'][-1][1:])
    
        f_map = dict(zip([x[0] for x in data['feature']], range(first.id, last.id + 1)))
        data['read_count'] = [
                [source_map[row[0]]] + [f_map[row[1]]] + row[2:] for row in data['read_count']
            ]
        t = ["read_source_id", "feature_id", "count", "span", "extend_up", "extend_down"]
    
        self.engine.execute(
            FeatureReadCount.__table__.insert(),
            [dict(zip(t, row)) for row in data['read_count']]
            )
    
        if data['evidence']:
            t = ["name","source"]
            result = self.engine.execute(
                Evidence.__table__.insert(),
                [dict(zip(t, row[1:])) for row in data['evidence']]
                )
    
            self.session.commit()
            first = self.fetch_evidence(data['evidence'][0][1:])
            last = self.fetch_evidence(data['evidence'][-1][1:])
    
            ev_map = dict(zip([x[0] for x in data['evidence']], range(first.id, last.id + 1)))
    
            data['feature_evidence'] = [
                    [f_map[row[0]], ev_map[row[1]]] for row in data['feature_evidence']
                ]
    
            t = ["feature_id", "evidence_id"]
            self.engine.execute(
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
                raise ValueError("Different chromosomes!")
            if e2[1] <= e1[2]:
                sys.stderr.write("{0} - {1}\n".format(e1, e2))
                raise ValueError("exons overlap, or in wrong order")
            if e1[3] != e2[3]:
                sys.stderr.write("{0} - {1}\n".format(e1, e2))
                raise ValueError("strands don't match")
       
        chrom = exons[0][0]
        strand = exons[0][-1]
        
        evidence = get_or_create(self.session, Evidence,
                name = name,
                source=source)

        seqs = []
        for exon in exons:
            seq = ""
            real_seq = ""
            if self.index:
                seq = ""
                try:                    
                    seq = self.index.get_sequence(chrom, exon[1] - 20, exon[2] + 20, strand)
                    real_seq = seq[20:-20]
                except Exception:
                    real_seq = self.index.get_sequence(chrom, exon[1], exon[2], strand)
                seqs.append(seq)
            
            exon = get_or_create(self.session, Feature,
                             chrom = chrom,
                             start = exon[1],
                             end = exon[2],
                             strand = strand,
                             ftype = "exon",
                             seq = real_seq
                             ) 
            exon.evidences.append(evidence)

        splice_donors = []
        splice_acceptors = []
        for i,(start,end) in enumerate([(e1[2], e2[1]) for e1, e2 in zip(exons[0:-1], exons[1:])]):
            self.logger.debug("%s %s %s %s", chrom, start, end, strand)
            sj = get_or_create(self.session, Feature,
                             chrom = chrom,
                             start = start,
                             end = end,
                             strand = strand,
                             ftype = "splice_junction"
                             )
            sj.evidences.append(evidence)
            
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
            self.logger.warning("Skipping %s, splicing not OK!", name)
            self.session.rollback()
        else:
            self.session.commit()
    
    def get_features(self, ftype=None, chrom=None, eager=False, 
            min_length=None, max_length=None, evidence=0):
        
        query = self.session.query(Feature)
        query = query.filter(Feature.flag.op("IS NOT")(True))
            
        # pre-fetch associated read counts
        if eager:
            query = query.options(subqueryload('read_counts'))
        
        if chrom:
            query = query.filter(Feature.chrom == chrom)
        if ftype:
            query = query.filter(Feature.ftype == ftype)
        features = [f for f in query]
        
        # length filters
        if max_length:
            features = [f for f in features if
                    f.length <= max_length or len(f.evidences) >= evidence]

        if min_length:
            features = [f for f in features if
                    f.length >= min_length or len(f.evidences) >= evidence]
        
        return features

    def get_exons(self, chrom=None, eager=False, min_length=None, 
            max_length=None, evidence=0):
        
        return self.get_features(ftype="exon", chrom=chrom, eager=eager, 
                min_length=min_length, max_length=max_length, evidence=evidence)


    def get_splice_junctions(self, chrom=None, ev_count=None, read_count=None, max_reads=None, eager=False):
                
        features = []
        if ev_count and read_count:
            # All splices with no read, but more than one evidence source
            fs = self.session.query(Feature).\
                    filter(Feature.flag.op("IS NOT")(True)).\
                    filter(Feature.ftype == "splice_junction").\
                    filter(Feature.chrom == chrom).\
                    outerjoin(FeatureReadCount).\
                    group_by(Feature).\
                    having(func.sum(FeatureReadCount.count) < read_count)

            for splice in fs:
                self.logger.debug("Considering %s", splice)        
                for evidence in splice.evidences:
                    self.logger.debug(str(evidence))
                if len(splice.evidences) >= ev_count:
                    features.append(splice)
                else:
                    self.logger.debug("not enough evidence for {}".format(splice))
            
            fs = self.session.query(Feature).\
                    filter(Feature.flag.op("IS NOT")(True)).\
                    filter(Feature.ftype == "splice_junction").\
                    filter(Feature.chrom == chrom).\
                    outerjoin(FeatureReadCount).\
                    group_by(Feature).\
                    having(func.sum(FeatureReadCount.count) == None)

            for splice in fs:
                self.logger.debug("Considering %s (no reads)", splice)        
                for evidence in splice.evidences:
                    self.logger.debug(str(evidence))
                if len(splice.evidences) >= ev_count:
                    features.append(splice)
                else:
                    self.logger.debug("not enough evidence for {}".format(splice))

            # All splices with more than x reads
            fs = self.session.query(Feature).\
                    filter(Feature.flag.op("IS NOT")(True)).\
                    filter(Feature.ftype == "splice_junction").\
                    filter(Feature.chrom == chrom).\
                    outerjoin(FeatureReadCount).\
                    group_by(Feature).\
                    having(func.sum(FeatureReadCount.count) >= read_count)
            for f in fs:
                self.logger.debug("Considering %s (reads)", f)        
                features.append(f)
            #features += [f for f in fs]
        
        elif max_reads:
            fs = self.session.query(Feature).\
                    filter(Feature.flag.op("IS NOT")(True)).\
                    filter(Feature.ftype == "splice_junction").\
                    filter(Feature.chrom == chrom).\
                    outerjoin(FeatureReadCount).\
                    group_by(Feature).\
                    having(func.sum(FeatureReadCount.count) == None)

            features += [f for f in fs if len(f.evidences) > 0]
            fs = self.session.query(Feature).\
                    filter(Feature.flag.op("IS NOT")(True)).\
                    filter(Feature.ftype == "splice_junction").\
                    filter(Feature.chrom == chrom).\
                    outerjoin(FeatureReadCount).\
                    group_by(Feature).\
                    having(func.sum(FeatureReadCount.count) < max_reads)
            features += [f for f in fs if len(f.evidences) > 0]  
        else:
            features = self.get_features(ftype="splice_junction", chrom=chrom, eager=eager)
        return features 

    def get_longest_3prime_exon(self, chrom, start5, strand):
        if strand == "+":
        
            q = self.session.query(Feature).\
                filter(Feature.flag.op("IS NOT")(True)).\
                filter(Feature.ftype == "exon").\
                filter(Feature.chrom == chrom).\
                filter(Feature.strand == strand).\
                filter(Feature.start == start5).\
                order_by(Feature.end)
            return q.all()[-1]
        else:
            q = self.session.query(Feature).\
                filter(Feature.ftype == "exon").\
                filter(Feature.chrom == chrom).\
                filter(Feature.strand == strand).\
                filter(Feature.end == start5).\
                order_by(Feature.end)
            return q.all()[0]

    
    def get_long_exons(self, chrom, l, evidence):
        query = self.session.query(Feature)
        query = query.filter(Feature.flag.op("IS NOT")(True))
        query = query.filter(Feature.ftype == 'exon')
        query = query.filter(Feature.chrom == chrom)
        query = query.filter(Feature.end - Feature.start >= l)
        return [e for e in query if len(e.evidences) <= evidence]
    
    def filter_repeats(self, chrom, rep):
        """ Flag all exons that overlap with a specified fraction
        with a repeat track
        """

        self.logger.warn("Filtering repeats: %s with fraction %s", 
                os.path.basename(rep["path"]), rep["fraction"])
        
        exons = self.get_features("exon", chrom)
        exon_tabix = exons_to_tabix_bed(exons) 
        
        overlap_it = tabix_overlap(exon_tabix, rep["tabix"], chrom, rep["fraction"])
        exon_ids = [int(iv[3]) for iv in overlap_it]
        
        chunk = 20
        for i in range(0, len(exon_ids), chunk):
            self.logger.warn("Filtering %s", exon_ids[i:i + chunk])
            self.session.query(Feature).\
                    filter(Feature.id.in_(exon_ids[i:i + chunk])).\
                    update({Feature.flag:True}, synchronize_session=False)
            self.session.commit()
            

        #fobj = TabixIteratorAsFile(tabixfile.fetch(chrom))
        #for line in fobj:
        #    print line

    def filter_evidence(self, chrom, source, experimental):
        self.logger.debug("Filtering %s", source)
        #query = self.session.query(Feature).\
        #        update({Feature.flag:False}, synchronize_session=False)
        #self.session.commit()

        # Select all features that are supported by other evidence
        n = self.session.query(Feature.id).\
                        join(FeatureEvidence).\
                        join(Evidence).\
                        filter(Evidence.source != source).\
                        filter(Evidence.source not in experimental).\
                        filter(Feature.chrom == chrom).\
                        subquery("n")
       
        # Select the total number of transcript from this source 
        # per feature
        s = self.session.query(Feature.id, func.count('*').label('total')).\
               join(FeatureEvidence).\
               join(Evidence).\
               filter(Evidence.source == source).\
               filter(Feature.chrom == chrom).\
               group_by(Feature.id).\
               subquery("s")

        # Select all features where support from this source
        # is only 1 transcript and which is not supported by any 
        # other sources
        a = self.session.query(Feature.id).filter(and_(
            Feature.id == s.c.id,
            s.c.total == 1)).\
            filter(Feature.id.notin_(n)).\
            subquery("a")

        #ids = [i[0] for i in query]
        
        # Flag features
        self.session.query(Feature).\
                filter(Feature.id.in_(a)).\
                update({Feature.flag:True}, synchronize_session=False)
        self.session.commit()
        
    def get_read_statistics(self, chrom, fnames, name, span="all", extend=(0,0), nreads=None):

        if span not in ["all", "start", "end"]:
            raise Exception("Incorrect span: {}".format(span))
        
        tmp = NamedTemporaryFile(delete=False)
        estore = {}
        self.logger.debug("Writing exons to file %s", tmp.name)
        exons =  self.get_exons(chrom)
        if len(exons) == 0:
            return
        
        for exon in exons:
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

            estr = "{}:{}-{}".format(exon.chrom, start, end)

            if estr in estore:
                estore[estr].append(exon)
            else:
                estore[estr] = [exon]
                tmp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
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
            self.logger.debug("Creating read_source for %s %s", name, fname)
            read_source = get_or_create(self.session, ReadSource, name=name, source=fname)
            self.session.commit() 
            #rmrepeats = False
            if fname.endswith("bam") and (not nreads or not nreads[i]):
                #rmrepeats = True
                self.logger.debug("Counting reads in %s", fname)
                read_source.nreads = read_statistics(fname)

            self.logger.debug("Getting overlap from %s", fname)
            result = get_binned_stats(tmp.name, fname, 1, rpkm=False, rmdup=False, rmrepeats=False)

            self.logger.debug("Reading results, save to exon stats")

            insert_vals = []
            for row in result:
                try:
                    vals = row.strip().split("\t")
                    e = "%s:%s-%s" % (vals[0], vals[1], vals[2])
                    c = float(vals[3])
                    for exon in estore[e]:
                        insert_vals.append([read_source.id, exon.id, c, span, extend[0], extend[1]])
                except:
                    self.logger.info("binned_stat line skipped: {}".format(row))
            t =  ["read_source_id", "feature_id", "count", "span", "extend_up", "extend_down"]
            result = self.engine.execute(
                    FeatureReadCount.__table__.insert(),
                    [dict(zip(t,row)) for row in insert_vals]
                    )
                
        tmp.close()

    def get_splice_statistics(self, chrom, fnames, name):
        if type("") == type(fnames):
            fnames = [fnames]

        for fname in fnames:
            self.logger.debug("Getting splicing data from %s", fname)
            read_source = get_or_create(self.session, 
                    ReadSource, name=name, source=fname)
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
            Feature.end == junction.start,
            Feature.ftype == "exon"
            )).\
            filter(Feature.flag.op("IS NOT")(True))

        
        right = self.session.query(Feature).filter(and_(
            Feature.chrom == junction.chrom,
            Feature.strand == junction.strand,
            Feature.start == junction.end,
            Feature.ftype == "exon"
            )).\
            filter(Feature.flag.op("IS NOT")(True))

        exon_pairs = []
        for e1 in left:
            for e2 in right:
                exon_pairs.append((e1, e2))
        return exon_pairs

    def clear_stats_cache(self):
        self.cache_feature_stats = {}
        self.cache_splice_stats = {}
    
    def feature_stats(self, feature, identifier):
        if "{}{}".format(feature, identifier) not in self.cache_feature_stats:
            q = self.session.query(FeatureReadCount, ReadSource).join(ReadSource)
            q = q.filter(FeatureReadCount.feature_id == feature.id)
            q = q.filter(ReadSource.name == identifier)
            self.cache_feature_stats["{}{}".format(feature, identifier)] = sum([x[0].count for x in q.all()])

        return self.cache_feature_stats["{}{}".format(feature, identifier)]

    def splice_stats(self, exon1, exon2, identifier):
        if "{}{}{}".format(self, exon1, exon2) not in self.cache_splice_stats:
            q = self.session.query(Feature)
            q = q.filter(Feature.ftype == "splice_junction")
            q = q.filter(Feature.chrom == exon1.chrom)
            q = q.filter(Feature.strand == exon1.strand)
            q = q.filter(Feature.start == exon1.end)
            q = q.filter(Feature.end == exon2.start)

            splice = q.first()
        
            self.cache_splice_stats["{}{}{}".format(self, exon1, exon2)] = self.feature_stats(splice, identifier)

        return self.cache_splice_stats["{}{}{}".format(self, exon1, exon2)]

    def intron_splice_stats(self, intron, identifier):
        if "{}{}".format(self, intron) not in self.cache_splice_stats:
            q = self.session.query(Feature)
            q = q.filter(Feature.ftype == "splice_junction")
            q = q.filter(Feature.chrom == intron.chrom)
            q = q.filter(Feature.strand == intron.strand)
            q = q.filter(Feature.start == intron.start)
            q = q.filter(Feature.end == intron.end)

            splice = q.first()
        
            self.cache_splice_stats["{}{}".format(self, intron)] = self.feature_stats(splice, identifier)

        return self.cache_splice_stats["{}{}".format(self, intron)]
    
    def nreads(self, identifier):
        q = self.session.query(ReadSource)
        q = q.filter(ReadSource.name == identifier)

        return sum([s.nreads for s in q.all()])
    
    def get_splice_count(self, e1, e2):
        counts = self.session.query(func.sum(FeatureReadCount.count)).\
                join(Feature).\
                filter(Feature.chrom == e1.chrom).\
                filter(Feature.start == e1.end).\
                filter(Feature.end == e2.start).\
                filter(Feature.ftype == "splice_junction").\
                group_by(Feature.id).all()
        return sum([int(x[0]) for x in counts])
    
    def fetch_feature(self, f):
        """ Feature as list """
        chrom, start, end, strand, ftype, seq = f

        feature = self.session.query(Feature).\
                filter(Feature.chrom == chrom).\
                filter(Feature.start == start).\
                filter(Feature.end == end).\
                filter(Feature.strand == strand).\
                filter(Feature.ftype == ftype)
                
        if seq:
            feature = feature.filter(Feature.seq == seq)
        result = feature.first()
        return result
    
    def fetch_evidence(self, f):
        """ Feature as list """
        name, source = f
        evidence = self.session.query(Evidence).\
            filter(Evidence.name == name).\
            filter(Evidence.source == source)

        result = evidence.first()
        return result

    def get_transcript_statistics(self, exons):
        stats = []
        for exon in exons:
            q = self.session.query(ReadSource, 
                    func.sum(FeatureReadCount.count)).\
                    join(FeatureReadCount).\
                    join(Feature).\
                    filter(Feature.chrom == exon[0]).\
                    filter(Feature.start == exon[1]).\
                    filter(Feature.end == exon[2]).\
                    filter(Feature.strand == exon[3]).\
                    group_by(ReadSource.name)
            stats.append(dict([(row[0].name, int(row[1])) for row in q.all()]))
        return stats

                   

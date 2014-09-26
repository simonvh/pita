import pytest

@pytest.fixture
def three_exons():
    exons = [
        ["chr1", 100, 200, "+"],
        ["chr1", 300, 400, "+"],
        ["chr2", 100, 200, "-"],
    ]
    return exons

@pytest.fixture
def three_transcripts():
    transcripts = [
        ["t1", "annotation",
            [
            ["chr1", 100, 200, "+"],
            ["chr1", 300, 400, "+"],
            ["chr1", 500, 600, "+"],
            ]
        ],
        ["t2", "annotation",
            [
            ["chr1", 50, 200, "+"],
            ["chr1", 300, 400, "+"],
            ]
        ],
        ["t3", "annotation",
            [
            ["chr1", 300, 400, "+"],
            ["chr1", 500, 800, "+"],
            ]
        ],
    ]
    return transcripts

@pytest.fixture
def two_transcripts():
    transcripts = [
        ["t1", "annotation",
            [
            ["chr1", 700, 900, "+"],
            ["chr1", 1100, 1200, "+"],
            ]
        ],
        ["t2", "annotation",
            [
            ["chr1", 1100, 1200, "+"],
            ["chr1", 1400, 1600, "+"],
            ]
        ],
    ]
    return transcripts

@pytest.fixture
def bam_file():
    return "tests/data/H3K4me3.bam"

@pytest.fixture
def transcripts():
    transcripts = [
            ["t1", "annotation",
                [
                    ["scaffold_1", 18070000, 18080000, "+"],
                    ["scaffold_1", 18200000, 18200100, "+"],
                    ["scaffold_1", 18250000, 18300000, "+"],
                ]
            ]
            
            ]
    return transcripts

@pytest.yield_fixture
def db(transcripts):
    from pita.annotationdb import AnnotationDb
    with AnnotationDb(new=True) as d:
        for name, source, exons in transcripts:
            d.add_transcript(name, source, exons)
        yield d
#scaffold_1  18070000    18080000    64
#scaffold_1  18080000    18200000    1092
#scaffold_1  18200000    18200100    1
#scaffold_1  18200100    18250000    318
#scaffold_1  18250000    18300000    300

def test_read_statistics(bam_file, db):
    db.get_read_statistics(bam_file, "test")
    exons = db.get_exons()
    counts = [e.read_counts[0].count for e in exons]
    assert [1,64,300] == sorted(counts)
    assert 5218 == db.nreads("test")

@pytest.fixture
def splice_file():
    return "tests/data/splice_data.bed"

def test_splice_statistics(db, splice_file):
    db.get_splice_statistics(splice_file, "test")
    splices = db.get_splice_junctions()
    counts = [s.read_counts[0].count for s in splices]
    assert 2 == len(splices)
    assert [4,20] == counts

def test_get_weight(db, bam_file, splice_file):
    db.get_read_statistics(bam_file, "H3K4me3")
    db.get_splice_statistics(splice_file, "RNAseq")
    from pita.dbcollection import DbCollection
    c = DbCollection(db)

    model = list(c.get_connected_models())[0][0]
    w = c.get_weight(model, "H3K4me3", "all")
    assert 365 == w
    w = c.get_weight(model, None, "length")
    assert 60100 == w
    w = c.get_weight(model, "H3K4me3", "rpkm")
    assert abs(1163.9 - w) < 0.1
    w = c.get_weight(model, "H3K4me3", "weighted")
    assert abs(0.01821963394342762 - w) < 0.0001
    w = c.get_weight(model, "H3K4me3", "total_rpkm")
    assert abs(4292.832 - w) < 0.1
    w = c.get_weight(model, "H3K4me3", "mean_exon")
    assert abs(1430.944 - w) < 0.1
    w = c.get_weight(model, "RNAseq", "splice")
    assert 24 == w
    w = c.get_weight(model, "H3K4me3", "first")
    assert 64 == w
    w = c.get_weight(model, None, "evidence")
    assert 1 == w

def test_get_junction_exons(db):
    splices = db.get_splice_junctions()
    splice = [s for s in splices if s.start == 18080000][0]

    exon_pairs = db.get_junction_exons(splice)
    assert 1 == len(exon_pairs)
    e1,e2 = exon_pairs[0]
    assert e1.start == 18070000
    assert e1.end == 18080000
    assert e2.start == 18200000
    assert e2.end == 18200100

def test_db_collection(db):    
    from pita.dbcollection import DbCollection
    c = DbCollection(db)

    for model in c.get_connected_models():
        print model

def test_get_long_exons(db):
    assert 0 == len(db.get_long_exons(100000))
    assert 1 == len(db.get_long_exons(50000))
    assert 2 == len(db.get_long_exons(10000))
    assert 3 == len(db.get_long_exons(50))


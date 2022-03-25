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
        [
            "t1",
            "annotation1",
            [
                ["chr1", 100, 200, "+"],
                ["chr1", 300, 400, "+"],
                ["chr1", 500, 600, "+"],
            ],
        ],
        [
            "t2",
            "annotation1",
            [
                ["chr1", 50, 200, "+"],
                ["chr1", 300, 400, "+"],
            ],
        ],
        [
            "t3",
            "annotation1",
            [
                ["chr1", 300, 400, "+"],
                ["chr1", 500, 800, "+"],
            ],
        ],
    ]
    return transcripts


@pytest.fixture
def two_transcripts():
    transcripts = [
        [
            "t1",
            "annotation2",
            [
                ["chr1", 700, 900, "+"],
                ["chr1", 1100, 1200, "+"],
            ],
        ],
        [
            "t2",
            "annotation2",
            [
                ["chr1", 1100, 1200, "+"],
                ["chr1", 1400, 1600, "+"],
            ],
        ],
    ]
    return transcripts


@pytest.fixture
def db(tmpdir):
    from pita.annotationdb import AnnotationDb

    conn = "sqlite:///{}/pita_test_database.db".format(tmpdir)
    db = AnnotationDb(conn=conn, new=True)
    return db


@pytest.fixture
def c(db):
    from pita.dbcollection import DbCollection

    c = DbCollection(db, [])
    return c


@pytest.fixture
def db_3t(db, three_transcripts):
    from pita.dbcollection import DbCollection

    for name, source, exons in three_transcripts:
        db.add_transcript(name, source, exons)
    c = DbCollection(db)
    return c


@pytest.fixture
def db_5t(db, two_transcripts, three_transcripts):
    from pita.dbcollection import DbCollection

    for name, source, exons in two_transcripts:
        db.add_transcript(name, source, exons)
    for name, source, exons in three_transcripts:
        db.add_transcript(name, source, exons)
    c = DbCollection(db, [])
    return c


def test_db_splice_junctions(db, two_transcripts, three_transcripts):
    from pita.dbcollection import DbCollection

    for name, source, exons in two_transcripts:
        db.add_transcript(name, source, exons, commit=True)
    for name, source, exons in three_transcripts:
        db.add_transcript(name, source, exons, commit=True)

    # No filter
    introns = db.get_splice_junctions(chrom="chr1")
    assert len(introns) == 4

    # Filter, but evidence 1 means get everything
    introns = db.get_splice_junctions(chrom="chr1", ev_count=1, read_count=1)
    assert len(introns) == 4

    # Strict filter, nothing as a result
    introns = db.get_splice_junctions(chrom="chr1", ev_count=3, read_count=10)
    assert len(introns) == 0

    # Filter, but keep annotation2
    introns = db.get_splice_junctions(
        chrom="chr1", ev_count=3, read_count=10, keep=["annotation2"]
    )
    assert len(introns) == 2

    # Filter, but keep annotation1
    introns = db.get_splice_junctions(
        chrom="chr1", ev_count=3, read_count=10, keep=["annotation1"]
    )
    assert len(introns) == 2

    # Filter, but keep annotation1
    introns = db.get_splice_junctions(
        chrom="chr1", ev_count=3, read_count=10, keep=["annotation1", "annotation2"]
    )
    assert len(introns) == 4


# def test_add_exon(db):
#    e = db.add_exon("chr1", 100, 200, "-")
#    assert "chr1:100-200" == str(e)
#    e = db.add_exon("chr1", 100, 200, "-")
#    assert "chr1:100-200" == str(e)


def test_add_transcripts(three_transcripts, db):
    for name, source, exons in three_transcripts:
        db.add_transcript(name, source, exons)

    assert 5 == len(db.get_exons())


# def test_get_initial_exons(db_3t):
#    assert 2 == len(db_3t.get_initial_exons())


def test_retrieve_models(db_5t):
    models = sorted(db_5t.get_best_variants([]), key=lambda x: len(x))
    assert 2 == len(models)
    assert 3 == len(models[0])
    assert 3 == len(models[1])


@pytest.fixture
def t1():
    return "tests/data/long_exons1.bed"


@pytest.fixture
def t2():
    return "tests/data/long_exons2.bed"


def test_long_exon_filter(db, t1, t2):
    from pita.dbcollection import DbCollection
    from pita.io import read_bed_transcripts

    for tname, source, exons in read_bed_transcripts(open(t1)):
        db.add_transcript("{0}{1}{2}".format("t1", "|", tname), source, exons)
    for tname, source, exons in read_bed_transcripts(open(t2)):
        db.add_transcript("{0}{1}{2}".format("t2", "|", tname), source, exons)

    c = DbCollection(db, [], chrom="chr1")
    c.filter_long(length=500, evidence=1)

    models = []
    for cluster in c.get_best_variants([]):
        models.append(cluster)

    assert [3, 5] == sorted([len(m) for m in models])


@pytest.fixture
def short_intron_track():
    return "tests/data/short_introns.bed"


# def test_short_intron_filter_merge(db, short_intron_track):
#    from pita.dbcollection import DbCollection
#    from pita.io import read_bed_transcripts
#
#    for tname, source, exons in read_bed_transcripts(open(short_intron_track)):
#        db.add_transcript("{0}{1}{2}".format("t1", "|", tname), source, exons)
#
#    c = DbCollection(db, chrom="chr1")
#    c.filter_short_introns(mode='merge')
#
#    exons = c.get_exons()
#    lens = sorted([len(e) for e in exons])
#    assert [100,100,100,110,1500,2000] == lens
#
# def test_short_intron_filter_delete(db, short_intron_track):
#    from pita.dbcollection import DbCollection
#    from pita.io import read_bed_transcripts
#
#    for tname, source, exons in read_bed_transcripts(open(short_intron_track)):
#        db.add_transcript("{0}{1}{2}".format("t1", "|", tname), source, exons)
#
#    c = DbCollection(db, chrom="chr1")
#
#    c.filter_short_introns(mode='delete')
#
#    exons = c.get_exons()
#    lens = sorted([len(e) for e in exons])
#    assert [100,100,100,1500,2000] == lens


@pytest.fixture
def variant_track():
    return "tests/data/many_paths.bed"


# def test_variants(db, variant_track):
#    from pita.dbcollection import DbCollection
#    from pita.io import read_bed_transcripts
#    from pita.util import model_to_bed
#
#    for tname, source, exons in read_bed_transcripts(open(variant_track)):
#         db.add_transcript("{0}{1}{2}".format("t1", "|", tname), source, exons)
#    c = DbCollection(db, [])
#
#    best_model = [m for m in  c.get_connected_models()][0][0]
#    cuts = [str(e) for e in c.get_node_cuts(best_model)]
#    assert ["chr1:800+900", "chr1:1400+1500"] == cuts
#
#    best_variant = c.get_best_variant(best_model, [{"weight":1,"type":"length","name":"length"}])
#    s = [str(s) for s in best_variant]
#    assert ["chr1:100+200", "chr1:400+700", "chr1:800+900", "chr1:1000+1300", "chr1:1400+1500", "chr1:1600+1900", "chr1:2000+2100"] == s
#
#

import pytest

@pytest.fixture
def bedfile():
    return "tests/data/overlapping_models.bed"

@pytest.fixture
def linkfile():
    return "tests/data/test_linkage.bed"

@pytest.fixture
def db():
    from pita.annotationdb import AnnotationDb
    db = AnnotationDb(conn="sqlite:///pita_test_database.db",
            new=True)
    return db

def test_get_overlapping_models(db, bedfile):
    from pita.io import read_bed_transcripts
    from pita.util import get_overlapping_models

    for tname, source, exons in read_bed_transcripts(open(bedfile), "test", 0):
        db.add_transcript("{0}{1}{2}".format("test", ":::", tname), source, exons)

    exons = db.get_exons("JGIv7b.000000226")
    
    assert 72 == len(get_overlapping_models(exons))

#def test_prune_cufflinks(linkfile):
#    from pita.collection import Collection
#    from pita.io import read_bed_transcripts
#    from pita.util import get_overlapping_models
#    mc = Collection()
#    
#    for tname, source, exons in read_bed_transcripts(open(linkfile), "test", 0):
#        mc.add_transcript("{0}{1}{2}".format("test", ":::", tname), source, exons)
#
#    mc.prune()
#    
#    assert 5 == len([x for x in mc.get_connected_models()])

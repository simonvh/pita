import pytest

@pytest.fixture
def bedfile():
    return "tests/data/overlapping_models.bed"

def test_get_overlapping_models(bedfile):
    from pita.collection import Collection
    from pita.io import read_bed_transcripts
    from pita.util import get_overlapping_models
    mc = Collection()

    for tname, source, exons in read_bed_transcripts(open(bedfile), "test", 0):
        mc.add_transcript("{0}{1}{2}".format("test", ":::", tname), source, exons)

    exons = mc.get_exons("JGIv7b.000000226")
    
    assert 72 == get_overlapping_models(exons)
     

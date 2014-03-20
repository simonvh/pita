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

def test_create_collection():
    from pita.collection import Collection
    c = Collection()

def test_add_exon():
    from pita.collection import Collection
    c = Collection()
    e = c.add_exon("chr1", 100, 200, "-")
    assert "chr1:100-200" == str(e)
    e = c.add_exon("chr1", 100, 200, "-")
    assert "chr1:100-200" == str(e)

def test_get_exons(three_exons):
    from pita.collection import Collection
    c = Collection()
    for e in three_exons:
        print e
        c.add_exon(*e)
    
    assert 3 == len(c.get_exons())
    assert 2 == len(c.get_exons("chr1"))
    assert 1 == len(c.get_exons("chr2"))

def test_add_transcripts(three_transcripts):
    from pita.collection import Collection
    c = Collection()
    for name, source, exons in three_transcripts:
        c.add_transcript(name, source, exons)

    assert 5 == len(c.get_exons())

def test_get_initial_exons(three_transcripts):
    from pita.collection import Collection
    c = Collection()
    for name, source, exons in three_transcripts:
        c.add_transcript(name, source, exons)

    assert 2 == len(c.get_initial_exons())
 
def test_retrieve_transcripts(three_transcripts):
    from pita.collection import Collection
    c = Collection()
    for name, source, exons in three_transcripts:
        c.add_transcript(name, source, exons)

    assert 4 == len(c.get_all_transcripts())
 
def test_retrieve_transcripts(three_transcripts, two_transcripts):
    from pita.collection import Collection
    c = Collection()
    for name, source, exons in three_transcripts + two_transcripts:
        c.add_transcript(name, source, exons)

    clusters = sorted(c.get_connected_models(), lambda x,y: cmp(len(x), len(y)))
    assert 2 == len(clusters)
    assert 1 == len(clusters[0])
    assert 4 == len(clusters[1])
    
@pytest.fixture
def t1():
    return "tests/data/long_exons1.bed" 

@pytest.fixture
def t2():
    return "tests/data/long_exons2.bed" 

def test_long_exon_filter(t1, t2):
    from pita.collection import Collection
    from pita.io import read_bed_transcripts
    from pita.util import model_to_bed

    c = Collection()

    for tname, source, exons in read_bed_transcripts(open(t1)):
        c.add_transcript("{0}{1}{2}".format("t1", "|", tname), source, exons)
    for tname, source, exons in read_bed_transcripts(open(t2)):
        c.add_transcript("{0}{1}{2}".format("t2", "|", tname), source, exons)
    
    c.filter_long()

    models = []
    for cluster in c.get_connected_models():
        for m in cluster:
            models.append(m)
    
    assert [1,3,5] == sorted([len(m) for m in models])

@pytest.fixture
def short_intron_track():
    return "tests/data/short_introns.bed" 

def test_long_exon_filter(short_intron_track):
    from pita.collection import Collection
    from pita.io import read_bed_transcripts
    from pita.util import model_to_bed

    c = Collection()

    for tname, source, exons in read_bed_transcripts(open(short_intron_track)):
        c.add_transcript("{0}{1}{2}".format("t1", "|", tname), source, exons)
    
    c.filter_short_introns(mode='merge')

    exons = c.get_exons()
    lens = sorted([len(e) for e in exons])
    assert [100,100,100,110,1500,2000] == lens
    
    c = Collection()

    for tname, source, exons in read_bed_transcripts(open(short_intron_track)):
        c.add_transcript("{0}{1}{2}".format("t1", "|", tname), source, exons)
    
    c.filter_short_introns(mode='delete')

    exons = c.get_exons()
    lens = sorted([len(e) for e in exons])
    assert [100,100,100,1500,2000] == lens

@pytest.fixture
def variant_track():
    return "tests/data/many_paths.bed" 

def test_variants(variant_track):
    from pita.collection import Collection
    from pita.io import read_bed_transcripts
    from pita.util import model_to_bed

    c = Collection()
    for tname, source, exons in read_bed_transcripts(open(variant_track)):
         c.add_transcript("{0}{1}{2}".format("t1", "|", tname), source, exons)

    best_model = [m for m in  c.get_connected_models()][0][0]
    cuts = [str(e) for e in c.get_node_cuts(best_model)]
    assert ["chr1:800-900", "chr1:1400-1500"] == cuts 
    
    best_variant = c.get_best_variant(best_model, [{"weight":1,"type":"length","name":"length"}])
    s = [str(s) for s in best_variant]
    assert ["chr1:100-200", "chr1:400-700", "chr1:800-900", "chr1:1000-1300", "chr1:1400-1500", "chr1:1600-1900", "chr1:2000-2100"] == s



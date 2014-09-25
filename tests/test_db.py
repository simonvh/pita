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

def test_create_db():
    from pita.db import AnnotationDb
    d = AnnotationDb()

def test_add_transcripts(three_transcripts):
    from pita.db import AnnotationDb
    d = AnnotationDb(new=True)
    for name, source, exons in three_transcripts:
        d.add_transcript(name, source, exons)

    assert 5 == len(d.get_exons())
    assert 5 == len(d.get_exons("chr1"))
    assert 0 == len(d.get_exons("chr2"))

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

#scaffold_1  18070000    18080000    64
#scaffold_1  18080000    18200000    1092
#scaffold_1  18200000    18200100    1
#scaffold_1  18200100    18250000    318
#scaffold_1  18250000    18300000    300

def test_read_statistics(bam_file, transcripts):
    from pita.db import AnnotationDb
    d = AnnotationDb(new=True)
    for name, source, exons in transcripts:
        d.add_transcript(name, source, exons)
   
    d.get_read_statistics(bam_file, "test")
    
    exons = d.get_exons()
    counts = [e.read_counts[0].count for e in exons]
    assert [1,64,300] == sorted(counts) 
    assert 1 == 0

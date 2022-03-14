import pytest

@pytest.fixture
def models_to_compare():
    from pita.io import read_bed_transcripts
    f1 = open("tests/data/annotation1.bed")
    f2 = open("tests/data/annotation2.bed")
    
    t1 = read_bed_transcripts(f1)
    t2 = read_bed_transcripts(f2)
    
    result = {}
    f = open("tests/data/ann1_vs_ann2.txt")
    f.readline() # header
    for line in f.readlines():
        vals = line.strip().split("\t")
        result.setdefault(vals[0], {})
        result[vals[0]][vals[1]] = [float(x) for x in vals[2:]]

    return t1, t2, result

def test_compare_annnotation(models_to_compare):
    from pita.compare import compare_annotation

    t1, t2, ref = models_to_compare

    result = compare_annotation(t1, t2)

    for gene1 in result.keys():
        for gene2 in result[gene1].keys():
            vals = result[gene1][gene2]
            ref_vals = ref[gene1][gene2]
            print(gene1, gene2, vals, ref_vals)
            for n, nref in zip(vals, ref_vals):
                assert abs(n - nref) < 0.0001

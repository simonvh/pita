import pytest

@pytest.fixture
def data():
    ret = []

    # No overlap
    ret.append([
        [0, 100, 200],
        [50, 50, 50],
        [0, 100, 200],
        [50, 50, 50],
    ])
    
    # Simple 2-exon overlap
    ret.append([
        [0, 100, 200],
        [101, 50, 50],
        [0, 200],
        [150, 50],
    ])
    # All exons overlap
    ret.append([
        [0, 100, 200],
        [101, 101, 50],
        [0],
        [250],
    ])
    # More than one overlap 
    ret.append([
        [0, 100, 200, 300, 400, 500, 600],
        [101, 101, 50, 50, 101, 101, 50],
        [0, 300, 400],
        [250, 50, 250],
    ])
    return ret

def test_retrieve_transcripts(data):
    from pita.io import merge_exons

    for starts, sizes, new_starts, new_sizes in data:
        assert new_starts, new_sizes == merge_exons(starts, sizes)
     

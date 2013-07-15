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
    return ret

def test_retrieve_transcripts(data):
    from pita.io import merge_exons

    for starts, sizes, new_starts, new_sizes in data:
        assert new_starts, new_sizes == merge_exons(starts, sizes)
     

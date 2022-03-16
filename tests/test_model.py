import pytest

from pita.annotationdb import AnnotationDb
from pita.model import load_chrom_data, get_chrom_models
from pita.util import model_to_bed


@pytest.fixture
def conn(tmpdir):
    conn = "sqlite:///{}/pita_test_database.db".format(tmpdir)
    return conn


def test_load_chrom_data(conn):
    data = []
    for ftype in ["bed", "gtf"]:
        data.append(
            [
                ftype,
                f"tests/data/gencode_test/test_gene.{ftype}",
                f"tests/data/gencode_test/test_gene.{ftype}.gz",
                ftype,
                2,
            ]
        )

    for dtype in data[:1], data[1:], data:
        db = AnnotationDb(new=True, conn=conn)

        load_chrom_data(
            conn,
            True,
            "chr9",
            dtype,
            [],
            "tests/data/genome/hg38_sample/hg38_sample.fa",
        )

        result = [
            gene
            for gene in get_chrom_models(
                conn, "chr9", [{"name": "orf", "type": "orf", "weight": 2}]
            )
        ]
        bed = model_to_bed(result[0][1])
        valid_result = "chr9\t14474\t29739\tchr9:14474-29739\t600\t-\t14806\t29662\t0,0,0\t11\t466,69,153,159,202,136,137,147,112,154,138,\t0,606,1434,2243,2490,2869,3244,3553,3906,10376,15127,"

        assert bed == valid_result

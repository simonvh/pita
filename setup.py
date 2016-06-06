from distutils.core import setup
from setuptools.command.test import test as TestCommand
import sys

VERSION = "1.76"
DESCRIPTION = """
pita - pita improves transcript annotation
"""
class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup (name = 'pita',
        version = VERSION,
        description = DESCRIPTION,
        author='Simon van Heeringen',
        author_email='s.vanheeringen@ncmls.ru.nl',
        license='MIT',
        url = 'https://github.com/simonvh/pita',
        download_url = 'https://github.com/simonvh/pita/tarball/'+VERSION,
        packages=[
            'pita'
        ],
        scripts=[
            "scripts/pita",
            "scripts/pita_utr",
            "scripts/bed12togff3",
            "scripts/gff3tobed12",
            "scripts/flatbread",
            "scripts/breadcrumb",
        ],
        data_files=[],
        tests_require=['pytest'],
      #  dependency_links = [
       #     'https://github.com/simonvh/gimmemotifs/archive/0.8.7.tar.gz#egg=gimmemotifs-0.8.7',
      #      ],
        install_requires=[
                        "numpy",
                        "SQLAlchemy",
                        "gimmemotifs > 0.8.6",
                        "pysam >= 0.9",
                        "pyyaml",
                        "HTSeq",
                        "bcbio-gff",
                        "biopython",
                        "networkx >= 1.10",
                        ],
        cmdclass = {'test': PyTest},
)

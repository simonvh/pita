pita improves transcript annotation 
===================================

Pipeline to improve transcript annotation based on RNA-seq and ChIP-seq data.

At the moment, it's still not in a stable state. It will be updated to annotate the Xenopus laevis genome based on experimental data.

Prerequisites
------------
The following Python modules are required:

* GFF parser - http://github.com/chapmanb/bcbb/tree/master/gff
* Biopython - http://biopython.org/
* pysam (>= 0.7.4)
* pyyaml
* networkx (>= 1.9)
* GimmeMotifs - http://github.com/simonvh/gimmemotifs
* HTSeq - http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
* numpy

Installation
------------

    # install prerequisites
    git clone git@bitbucket.org:simonvh/pita.git
    cd pita
    python setup.py test
    sudo python setup.py install
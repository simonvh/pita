pita improves transcript annotation 
===================================

Pipeline to improve transcript annotation based on RNA-seq and ChIP-seq data.

At the moment, it's still a mess. It will be updated to annotate the Xenopus laevis genome with based on experimental data.

Prerequisites
------------
The following Python modules are required:

* GFF parser - http://github.com/chapmanb/bcbb/tree/master/gff
* Biopython - http://biopython.org/
* pysam (>= 0.7.4)
* pyyaml
* networkx (>= 1.9)

Installation
------------

    # install prerequisites
    git clone git@bitbucket.org:simonvh/pita.git
    cd pita
    sudo python setup.py install

pita improves transcript annotation 
===================================

[![Build Status](https://travis-ci.org/simonvh/pita.svg?branch=master)](https://travis-ci.org/simonvh/pita)
[![Code Health](https://landscape.io/github/simonvh/pita/master/landscape.svg?style=flat)](https://landscape.io/github/simonvh/pita/master)


Pipeline to improve transcript annotation based on RNA-seq and ChIP-seq data.

The current version has been used to annotate the Xenopus laevis genome based on experimental data.

However, it is not yet easy to install and use as the documentation is incomplete. 
In addition  the tools have not been thoroughly tested on a clean installation, 
which means I'm not sure all dependencies have been correctly specified.

Prerequisites
------------
The following Python modules are required:

* GFF parser - http://github.com/chapmanb/bcbb/tree/master/gff
* Biopython - http://biopython.org/
* pysam ( >= 0.9)
* pyyaml
* networkx (== 1.9)
* GimmeMotifs - http://github.com/simonvh/gimmemotifs
* HTSeq - http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
* numpy

Installation
------------

    # install prerequisites
    git clone git@github.com/simonvh/pita.git
    cd pita
    python setup.py test
    sudo python setup.py install


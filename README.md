# pita improves transcript annotation 

[![Build Status](https://travis-ci.org/simonvh/pita.svg?branch=master)](https://travis-ci.org/simonvh/pita)
[![Code Health](https://landscape.io/github/simonvh/pita/master/landscape.svg?style=flat)](https://landscape.io/github/simonvh/pita/master)


Pipeline to improve transcript annotation based on RNA-seq and ChIP-seq data.

It has been used to annotate the Xenopus laevis genome based on experimental data.

> Session et al. Genome evolution in the allotetraploid frog Xenopus laevis. Nature. 2016 Oct 20;538(7625):336-343. [doi: 10.1038/nature19840](http://dx.doi.org/10.1038/nature19840).

However, it is not yet easy to use as the documentation is incomplete. 

## Installation

Clone this repository:

```
$ git clone git@github.com/simonvh/pita.git
$ cd pita
```

Install dependencies via [conda](https://conda.io/):

```
$ conda env create -f environment.yml
```

Activate the environment:

```
$ source activate pita
```

Test & install:

```
$ python setup.py test
$ python setup.py install
```

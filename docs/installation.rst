Installation
============

Pita runs on Linux. Definitely not on Windows, sorry. Mac OS X
should work in theory, but as I don’t have the means to test this, I’m
not completely sure.


Required packages
-----------------


- Python v2.7 (Python 3 is not supported)

- Pandas v9.0

  Python package, providing data structures and analysis tools (installable with pip).

- Fluff v1.62

  Fluff is a python package that contains several scripts to produce publication-quality figures for next g$
  The source is available on github `<https://github.com/simonvh/fluff>`_


Using pip
---------

The most straightforward option is installing with pip. Pita is not yet officially hosted on PyPi (this option will be available soon) but can be installed using pip with the following command:

::

	$ sudo pip install git+https://github.com/simonvh/pita

Using pip in a virtualenv
--------------------------
If you don't have root acces (for example when working on a server) pita can be easily installed in a virtual environment using the following commands.

::

	$ virtualenv pitaEnv
	$ source pitaEnv/bin/activate
	$ pip install git+https://github.com/simonvh/pita


Installation from source
------------------------

::

	# install prerequisites
	$ git clone git@github.com/simonvh/pita.git
	$ cd pita
	$ python setup.py test
	$ sudo python setup.py install




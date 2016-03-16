Scripts
=======


bam2splicecount
---------------

bam2splicecount takes a sorted bam file as input, extracts the splice sites 
(e.g. exon boundaries) and counts these sites for multiple occurences in the following format.
The script can consume a lot of memory dependent on the size of the BAM file.

  =====  =======  =======   =======
   Chr    Start    Stop      Count
  =====  =======  =======   =======
  Chr01   23000    23001       2
  =====  =======  =======   =======

**Usage:**

::

	$ bam2splicecount <bamfile>

bam2splicecount takes only one argument: <bamfile>.
The output is redirected to stdout.


bed12togff3
-----------

The name explains it all, this script converts bed12 files to the gff3 format.

**Usage:**

::

	$ bed12togff3 <bed12file>

bed12togff3 takes only one argument: <bed12file>.
The output is redirected to stdout.

breadcrumb
----------

Coming soon..


flatbread
---------

Flatbread is a script for creating test datasets. The script accepts a yaml 
configuration file as an argument, in the same way as Pita does, along with a bed 
file containing regions. Flatbread then processes al the data and annotation files
present in the yaml configuration and filters out only the regions specified in the
bed file and writes them to new files. The result is data bam files and annotation files
for only the regions present in the bed file, that are much smaller and can be 
used for testing of the pipeline. The second result is a new yaml file with references
to the new files that can be directly used in pita.

**Usage:**

::

	$ flatbread [-h] [-c CONFIGFILE] [-b BEDFILE] [-o OUTPUT]

- ``-c`` Input configuration file in yaml format. Containing references to all the
data and annotation files that have to be filtered.

- ``-b`` BED file with genomic regions. These are the regions that will be present
in the output of this script.

- ``-o`` Output name. A directory with this name will be created containing all the
output files.


gff3tobed12
-----------

The name explains it all, this script converts gff3 files to the bed12 format.

**Usage:**

::

	$ gff3tobed12 <gff3>

The output is directed to stdout.

yamlMaker
---------

yamlMaker is a script to easily create a yaml file which can be used a input
for pita/flatbread/breadcrumb. When ran, the script asks the user questions about
the data or annotation entries. The result is a configuration file in yaml format.

**Usage:**

::

	$ yamlMaker <output.yaml>  

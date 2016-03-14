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

bam2splicecount takes only one argument: <bamfile>.
The output is redirected to stdout.

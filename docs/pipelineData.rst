Pipeline data
============

Pita is capable of using various types of information. This section describes which data types are currently supported.


Raw data 
--------
Two types of data generally should be present when running the pipeline:

- RNA-seq
  
  Mapped RNA-seq data in bam files.

- ChIP-seq

  Mapped ChIP-seq data in bam files. There are multiple types of ChIP-seq data. Pita was tested on H3K4 tri-methylations.


Annotation
----------

- **Splice sites**

  Splice sites counts for each bam file in the Raw data section, according to the following format

  =====  =======  =======   =======
   Chr    Start    Stop      Count
  =====  =======  =======   =======
  Chr01   23000    23001       2
  =====  =======  =======   =======
  
  Example script: `<https://github.com/simonvh/pita/blob/master/scripts/bam2splicecount>`_  

- **Assembled transcripts**
 
  Transcript assembled based on RNA-seq data. For testing we used the tool stringtie on each seperate RNA-seq bam file (guided by a reference file)
  and merged them after.

  ::

	$ stringtie input.bam -G /path/to/geneModels.gff3 -o output.gtf (on each bam file)
	$ stringtie --merge -G /path/to/geneModels.gff3 [list of gtf files]

- **ESTs & cDNA**
  
  These two data types can be used as an input. For testing we retrieved the sequences from sources such as xenbase and ensembl. 
  However just the sequences are not enough, we need the locations on the genome which requires mapping to the reference genome.
  We reccomend gmap to do this.

- **Protein sequences**
 
  Protein sequences can be very informative to predict coding regions. Pita is not able to directly use the amino acid sequences so some preprocessing is required.
  Using blat the amino acids can be mapped to the genome resulting in the location of the proteins. Blat does a good job in mapping, however the splice sites are
  usually not great. It is recommended to correct the splice sites with scipio. Scipio uses the result psl file from Blat as an input file to generate a 
  gff file with corrected splice sites. The resulting gff file can be used as an input file for Pita .

- **Existing gene models**

  Existing gene models can be used to guide the annotation process and usually do not need any preprocessing. 
  The prefered format is bed, but gff is supported as well. 



References
__________


- Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT & Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads Nature Biotechnology 2015, doi:10.1038/nbt.3122

- Wu, T. D., & Watanabe, C. K. (2005). GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics (Oxford, England), 21(9), 1859–75. http://doi.org/10.1093/bioinformatics/bti310

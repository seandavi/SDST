**************************************
Using the ``seqtool`` script
**************************************

Overview
========

The ``seqtools`` script is the entry point for command-line-oriented interaction with seqtools. 
Usage is pretty straightforward.  To get top-level help, simply type:

.. code-block:: bash
   
   seqtool -h



Usage Examples
==============

Demultiplexing FastQ Files
--------------------------

.. code-block:: bash

   $ seqtool utils demultiplex -h

This will return the help for the demultiplexer.

.. code-block:: bash

   usage: seqtool utils demultiplex [-h] [-1 READFILE1] [-2 READFILE2]
				    [-i INDEXFILE] [-j INDEXFILE2] [-x INDEX]

   optional arguments:
     -h, --help            show this help message and exit
     -1 READFILE1, --readFile1 READFILE1
			   read1 filename
     -2 READFILE2, --readFile2 READFILE2
			   read2 filename
     -i INDEXFILE, --indexFile INDEXFILE
			   index Filename
     -j INDEXFILE2, --indexFile2 INDEXFILE2
			   index Filename
     -x INDEX, --index INDEX
			   The indexes, one per index or comma-separated pairs
			   for dual-barcode indexing



Fixing VarScan VCF Files
------------------------

.. code-block:: bash

   $ seqtool varscan fixVcf -h
   usage: seqtool varscan fixVcf [-h] [-f VARSCAN]

   optional arguments:
     -h, --help            show this help message and exit
     -f VARSCAN, --varscan VARSCAN
                           Filename of varscan VCF file [default=stdin]

Melt a VCF File to Tab-Delimited Text
-------------------------------------

A seqtool vcf subcommand, melt, can take a VCF file (including multisample VCFs) and output a
tab-delimited text file, expanding all the various INFO and FORMAT columns as well as the
snpEff output.  For snpEff, note that the highest impact variant is chosen (based on the 
snpEff classification scheme.   

.. code-block:: bash

   $ seqtool vcf melt -f PATEEM_DNA.snp_fixed.final.vcf.gz -o PATEEM_DNA.snp_fixed.final.txt



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
                                    [-i INDEXFILE] [-x INDEX]

   optional arguments:
     -h, --help            show this help message and exit
     -1 READFILE1, --readFile1 READFILE1
                           read1 filename
     -2 READFILE2, --readFile2 READFILE2
                           read2 filename
     -i INDEXFILE, --indexFile INDEXFILE
                           index Filename
     -x INDEX, --index INDEX
                           The indexes, one per index


Fixing VarScan VCF Files
------------------------

.. code-block:: bash

   $ seqtools varscan fixVcf -h
   usage: seqtool varscan fixVcf [-h] [-f VARSCAN]

   optional arguments:
     -h, --help            show this help message and exit
     -f VARSCAN, --varscan VARSCAN
                           Filename of varscan VCF file [default=stdin]

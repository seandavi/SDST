======================
Working with VCF files
======================

Annotating VCF Files
====================

snpEff
------

Download data
^^^^^^^^^^^^^

.. code-block:: bash 

    java -jar /data/CCRBioinfo/biowulf/local/snpEff_3_0/snpEff.jar download -c /data/CCRBioinfo/biowulf/local/snpEff_3_0/snpEff.config GRCh37.66

Run effect prediction
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash 

    java -jar /data/CCRBioinfo/biowulf/local/snpEff_3_0/snpEff.jar eff -c /data/CCRBioinfo/biowulf/local/snpEff_3_0/snpEff.config GRCh37.66


snpSift
-------

Annotate with dbSNP
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash
 
    java -Xmx8g -jar /data/CCRBioinfo/biowulf/local/SnpSift_latest.jar annotate  /data/CCRBioinfo/public/GATK/bundle/1.5/hg19/dbsnp_135.hg19.excluding_sites_after_129.vcf tmp2.vcf > tmp2.dbsnp.vcf



Annotate with dbNSFP
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash 

    # Donload and uncompress database (you need to do this only once):
    # WARNING: The database is 3Gb when compressed and 30Gb uncompressed.
    wget http://sourceforge.net/projects/snpeff/files/databases/dbNSFP2.0b3.txt.gz
    gunzip dbNSFP2.0b3.txt.gz

    java -jar SnpSift.jar dbnsfp /data/CCRBioinfo/public/snpEff/data/dbNSFP2.0b3.txt myFile.vcf > myFile.annotated.vcf


Filtering
^^^^^^^^^

The snpSift package allows very flexible filtering options.  An example is given here:

.. code-block:: bash 

    java -Xmx8g -jar /data/CCRBioinfo/biowulf/local/SnpSift_latest.jar filter '(na ID) & (ID =~ 'COSM') & !( ID =~ 'rs')' -f 

Melting to Tab-delimited Text
-----------------------------

VCF files are quite difficult to read and filter in something like Excel.  The term, "vcf melting", refers to taking a VCF file and pulling out the various parts into separate tab-delimited columns.  See ``seqtool vcf melt`` for details of command-line operation or seqtools.vcf.vcfMelt for API usage.




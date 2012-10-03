Structural Variation
====================

GASV on Helix/Biowulf
---------------------

.. code-block:: bash

   /usr/local/bin/java64 -Xmx8g -jar /data/CCRBioinfo/biowulf/local/src/gasv/bin/BAMToGASV.jar 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam -MAPPING_QUALITY 30 -CUTOFF_LMINLMAX SD=3 -LIBRARY_SEPARATED all


Output from this step looks like this:


.. code-block:: bash

   -rw-r----- 1 sdavis2 sdavis2  39G Jul 27 15:54 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam
   -rw-r--r-- 1 sdavis2 sdavis2 165M Jul 30 12:39 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam_all.deletion
   -rw-r--r-- 1 sdavis2 sdavis2 1.1M Jul 30 12:39 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam_all.inversion
   -rw-r--r-- 1 sdavis2 sdavis2 192M Jul 30 12:39 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam_all.divergent
   -rw-r--r-- 1 sdavis2 sdavis2 2.2M Jul 30 12:39 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam_all.translocation
   -rw-r--r-- 1 sdavis2 sdavis2    0 Jul 30 12:39 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam_all.insertion
   -rw-r--r-- 1 sdavis2 sdavis2   32 Jul 30 12:39 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam.info
   -rw-r--r-- 1 sdavis2 sdavis2  331 Jul 30 12:39 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam.gasv.in
   

The next step is to use GASV to cluster and call the SVs:

.. code-block:: bash

   /usr/local/bin/java64 -Xmx20g -jar /data/CCRBioinfo/biowulf/local/src/gasv/bin/GASV.jar --fast --batch 5_BD13EMACXX.361_BUSTARD-2012-07-05.fq.hg19.bwa.sorted.bam.gasv.in

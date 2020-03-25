PPanGGOLiN : Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors
========================================================================================================

PPanGGOLiN is a software suite used to create and manipulate prokaryotic pangenomes from a set of either genomic DNA sequences or provided genome annotations. It is designed to scale up to tens of thousands of genomes. It has the specificity to partition the pangenome using a statistical approach rather than using fixed thresholds which gives it the ability to work with low-quality data such as Metagenomic Assembled Genomes (MAGs) or Single-cell Amplified Genomes (SAGs) thus taking advantage of large scale environmental studies and letting users study the pangenome of uncultivable species.

PPanGGOLiN builds pangenomes through a graphical model and a statistical method to partition gene families in persistent, shell and cloud genomes. It integrates both information on protein-coding genes and their genomic neighborhood to build a graph of gene families where each node is a gene family and each edge is a relation of genetic contiguity. The partitioning method promotes that two gene families that are consistent neighbors in the graph are more likely to belong to the same partition. It results in a Partitioned Pangenome Graph (PPG) made of persistent, shell and cloud nodes drawing genomes on rails like a subway map to help biologists navigate the great diversity of microbial life.

|qual| |cov| |bioconda| |plat| |version|

.. |qual| image:: https://api.codacy.com/project/badge/Grade/a24bff9354504a3294f4acf70681765a
.. |cov| image:: https://api.codacy.com/project/badge/Coverage/806fdcd8d04a469e8233728780576160
.. |plat| image:: https://anaconda.org/bioconda/ppanggolin/badges/platforms.svg
   :target: https://anaconda.org/bioconda/ppanggolin
.. |version| image:: https://anaconda.org/bioconda/ppanggolin/badges/version.svg
   :target: https://anaconda.org/bioconda/ppanggolin
.. |bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target: http://bioconda.github.io/recipes/ppanggolin/README.html


.. image:: images/logo.png
    :align: center


Installation
============


PPanGGOLiN is easily installed via conda
You will need the following conda channels if you don't have them already:

.. code:: bash

	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda config --add channels r

And then you can just run :

.. code:: bash

	conda install -c bioconda ppanggolin

Quick usage
===========

PPanGGOLiN has a minimal command for the non-expert users :

.. code:: bash

	ppanggolin workflow --fasta ORGANISMS_FASTA_LIST

It uses parameters that we found to be generally the best when working with species pangenomes.

The file ORGANISMS_FASTA_LIST is a tsv-separated file with the following organisation :
	1. The first column contains a unique organism name **(without whitespace)**
	2. The second column the path to the associated FASTA file
	3. Circular contig identifiers are indicated in the following columns
	4. Each line represents an organism

An `example <https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.fasta.list>`_ with 50 *Chlamydia trachomatis* genomes can be found in the testingDataset/ directory.

You can also give PPanGGOLiN your own annotations using .gff or .gbff/.gbk files instead of .fasta files, such as the ones provided by prokka using the following command :

.. code:: bash

	ppanggolin workflow --anno ORGANISMS_ANNOTATION_LIST

Another `example <https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list>`_ of such a file can be found in the testingDataset/ directory.

Both of those commands write several output files and graphics. Most notably a HDF-5 (pangenome.h5) file is written. It can be used as input for any of the subcommands to rerun parts of the analysis with different parameters, write and draw different representations of the pangenome or run additional analysis with PPanGGOLiN.

A minimum of 5 genomes is generally required to perform a pangenomics analysis using the traditional *core genome*/*accessory genome* paradigm. It is advised to use at least 15 genomes having genomic variations (and not only SNPs) to obtain robust results with the PPanGGOLiN statistical approach.

If you want to use personalized parameters for each subcommand most options should be self descriptive. If you want to know more about what each output file is, or briefly how each subcommand works you can check the `github wiki <https://github.com/labgem/PPanGGOLiN/wiki>`_

Issues, Questions, Remarks
==========================

If you have any question or issue with installing, using or understanding PPanGGOLiN, please do not hesitate to post an issue ! We cannot correct bugs if we do not know about them, and will try to help you the best we can.


Citation
========
If you use this tool for your research please cite :

Gautreau G et al. (2020) PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph. PLOS Computational Biology 16(3): e1007732. https://doi.org/10.1371/journal.pcbi.1007732

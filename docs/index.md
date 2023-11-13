% PPanGGOLiN documentation master file, created by
% sphinx-quickstart on Tue Sep 12 10:29:06 2023.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

# PPanGGOLiN: Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyrodigal/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/labgem/ppanggolin/actions)
[![License](https://anaconda.org/bioconda/ppanggolin/badges/license.svg)](http://www.cecill.info/licences.fr.html)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/ppanggolin?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/ppanggolin)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/labgem/ppanggolin/)
[![GitHub issues](https://img.shields.io/github/issues/labgem/ppanggolin.svg?style=flat-square&maxAge=600)](https://github.com/labgem/ppanggolin/issues)
[![Docs](https://img.shields.io/readthedocs/ppanggolin/latest?style=flat-square&maxAge=600)](https://ppanggolin.readthedocs.io)
[![Downloads](https://anaconda.org/bioconda/ppanggolin/badges/downloads.svg)](https://bioconda.github.io/recipes/ppanggolin/README.html#download-stats)

**PPanGGOLiN**
([Gautreau et al. 2020](https://doi.org/10.1371/journal.pcbi.1007732)) is a software suite used to create and manipulate prokaryotic pangenomes from a set of either genomic DNA sequences or provided genome annotations.
It is designed to scale up to tens of thousands of genomes.
It has the specificity to partition the pangenome using a statistical approach rather than using fixed thresholds which gives it the ability to work with low-quality data such as *Metagenomic Assembled Genomes (MAGs)* or *Single-cell Amplified Genomes (SAGs)* thus taking advantage of large scale environmental studies and letting users study the pangenome of uncultivable species.

**PPanGGOLiN** builds pangenomes through a graphical model and a statistical method to partition gene families in persistent, shell and cloud genomes.
It integrates both information on protein-coding genes and their genomic neighborhood to build a graph of gene families where each node is a gene family, and each edge is a relation of genetic contiguity.
The partitioning method promotes that two gene families that are consistent neighbors in the graph are more likely to belong to the same partition.
It results in a Partitioned Pangenome Graph (PPG) made of persistent, shell and cloud nodes drawing genomes on rails like a subway map to help biologists navigate the great diversity of microbial life.


Moreover, the panRGP method ([Bazin et al. 2020](https://doi.org/10.1093/bioinformatics/btaa792)) included in **PPanGGOLiN** predicts, for each genome, Regions of Genome Plasticity (RGPs) that are clusters of genes made of shell and cloud genomes in the pangenome graph.
Most of them arise from Horizontal gene transfer (HGT) and correspond to Genomic Islands (GIs). 
RGPs from different genomes are next grouped in spots of insertion based on their conserved flanking persistent genes.


Those RGPs can be further divided in conserved modules by panModule ([Bazin et al. 2021](https://doi.org/10.1101/2021.12.06.471380)). Those conserved modules correspond to groups of cooccurring and colocalized genes that are gained or lost together in the variable regions of the pangenome.

```{image} _static/logo.png
:alt: ppangolin logo
:align: center
:heigth: 300
:width: 300
```

# Installation

**PPanGGOLiN** is easily installed via conda. 
You will need the following conda channels if you don't have them already:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Then, you can just run :

```bash
conda install -c bioconda ppanggolin
```

# Quick usage

**PPanGGOLiN** integrates some workflows to build and analyse easily and rapidly a pangenome. 
These commands can be tuned with some parameters but are mostly automatic.
All workflow parameters are described [here](#step-by-step-section).

## Pangenome graph construction and partition

To build and partition a pangenome, you can use the following command:
```bash
ppanggolin workflow --fasta ORGANISMS_FASTA_LIST
```

It uses parameters that we found to be generally the best when working with species pangenomes.

The file ORGANISMS_FASTA_LIST is a tsv-separated file with the following organization :
1. The first column contains a unique organism name **(without space)**
2. The second column the path to the associated FASTA file
3. Circular contig identifiers are indicated in the following columns
4. Each line represents an organism

An [example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.fasta.list) with 50 *Chlamydia trachomatis* genomes can be found in the testingDataset/ directory.


You can also give **PPanGGOLiN** your own annotations using *.gff* or *.gbff/.gbk* files instead of *.fasta* files,
such as the ones provided by prokka using the following command :

```bash
ppanggolin workflow --anno ORGANISMS_ANNOTATION_LIST
```

Another [example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list) of such a file can be found in the testingDataset/ directory.

Both of those commands write several output files and graphics (more information [here](#output)). Most notably, an HDF-5 (pangenome.h5) file is written.
It can be used as input for any of the subcommands to rerun parts of the analysis with different parameters,
write and draw different representations of the pangenome or run additional analysis with **PPanGGOLiN**.

A minimum of 5 genomes is generally required to perform a pangenomics analysis using the traditional *core genome*/*accessory genome* paradigm.
It is advised to use at least 15 genomes having genomic variations (and not only SNPs) to obtain robust results with the **PPanGGOLiN** statistical approach.

If you want to use personalized parameters for each subcommand, most options should be self-descriptive.
If you want to know more about what each output file is, or briefly how each subcommand works,
you can check the [steb by step documentation](https://github.com/labgem/PPanGGOLiN/wiki)


## Region of plasticity detection

Furthermore, you can also predict genomic islands and cluster them into spots of insertion using the **panRGP** pipeline.
The usage is identical to the previous 'workflow' command:

```bash
ppanggolin panrgp --fasta ORGANISMS_FASTA_LIST
```

It will run more analyses after the pangenome has been partitioned. Further details are available [here](#panrgp) and in the [panRPG publication](https://doi.org/10.1093/bioinformatics/btaa792)

## Conserved module prediction
To detect the conserved modules in your pangenome, you can use the panModule workflow, as such:

```bash
ppanggolin panmodule --fasta ORGANISMS_FASTA_LIST
```

Further details can be found [here](#panmodule) and in the [panModule publication](https://doi.org/10.1101/2021.12.06.471380)


Alternatively, to run all the possible analysis that **PPanGGOLiN** can run, you can use:

```bash
ppanggolin all --fasta ORGANISMS_FASTA_LIST
```

Overall, ppanggolin has a lot of subcommands and possibilities.
Don't hesitate to check the command line help, and the [GitHub wiki](https://github.com/labgem/PPanGGOLiN/wiki) to see all the possible analysis, if you are missing a file you're looking for, or do not understand an output.
You can also raise an `issue` if you wish!

# Issues, Questions, Remarks
If you have any question or issue with installing,
using or understanding **PPanGGOLiN**, please do not hesitate to post an issue!
We cannot correct bugs if we do not know about them, and will try to help you the best we can.

# Citation
If you use this tool for your research, please cite:

Gautreau G et al. (2020) **PPanGGOLiN**: Depicting microbial diversity via a partitioned pangenome graph.
PLOS Computational Biology 16(3): e1007732. <https://doi.org/10.1371/journal.pcbi.1007732>

If you use this tool to study genomic islands, please cite:

Bazin et al., panRGP: a pangenome-based method to predict genomic islands and explore their diversity, Bioinformatics, Volume 36, Issue Supplement_2, December 2020, Pages i651–i658, <https://doi.org/10.1093/bioinformatics/btaa792>

If you use this tool to study modules, please cite:

Bazin et al., panModule: detecting conserved modules in the variable regions of a pangenome graph. biorxiv. <https://doi.org/10.1101/2021.12.06.471380>


```{toctree}
:caption: 'User Guide:'
:maxdepth: 2

user/Introduction
user/Installation
user/Basic-usage-and-practical-information
user/step-by-step-pangenome-analysis
user/Regions-of-Genome-Plasticity
user/Conserved-modules
user/Align
user/Genomic-context
user/projection
user/metadata
user/Outputs
```

```{toctree}
:caption: 'Developper Guide:'
:maxdepth: 2

dev/devRules
dev/git
dev/unitTest
dev/workflows
dev/buildDoc
```

# Indices and tables
[//]: # (- {ref}`ppanggolin package`)

- {ref}`genindex`

- {ref}`modindex`

- {ref}`search`

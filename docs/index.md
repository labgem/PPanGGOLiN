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


```{image} _static/logo.png
:alt: ppangolin logo
:align: center
:width: 437
```

<br>


**PPanGGOLiN**
([Gautreau et al. 2020](https://doi.org/10.1371/journal.pcbi.1007732)) is a software suite used to create and manipulate prokaryotic pangenomes from a set of either assembled genomic DNA sequences or provided genome annotations.
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


<!-- 
[//]: # (```{toctree})

[//]: # (:caption: 'Tutorial:')

[//]: # (:maxdepth: 1)

[//]: # ()
[//]: # (tutorial/prepEnv)

[//]: # (tutorial/inputData)

[//]: # (tutorial/workflows)

[//]: # (tutorial/AnalyseRGP)

[//]: # (```) -->

```{toctree}
:caption: 'User Guide:'
:maxdepth: 2

user/install
user/QuickUsage/quickAnalyses
user/practicalInformation
user/PangenomeAnalyses/pangenomeAnalyses
user/RGP/rgpAnalyses
user/Modules/moduleAnalyses
user/writeGenomes
user/align
user/projection
user/genomicContext
user/MSA
user/metadata
```

```{toctree}
:caption: 'Developper Guide:'
:maxdepth: 1

dev/contribute
dev/buildDoc
api/modules
```


# Indices and tables
[//]: # (- {ref}`ppanggolin package`)

- {ref}`genindex`

- {ref}`modindex`

- {ref}`search`

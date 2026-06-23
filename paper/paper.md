---
title: 'PPanGGOLiN V2: technical enhancement and extended functionalities for prokaryotic pangenome analysis'

tags:
- Python
- Bioinformatics
- Pangenomics
- Microbial genomics
- Mobile genetic elements
- Comparative genomics

authors:
- name: Jérôme Arnoux
  orcid: 0000-0003-2769-3006
  affiliation: 1
  equal-contrib: true
- name: Jean Mainguy
  orcid: 0009-0006-9160-9744
  affiliation: 1
  equal-contrib: true
- name: Adelme Bazin
  orcid: 0000-0002-5656-4708
  affiliation: 2
- name: Guillaume Gautreau
  orcid: 0000-0002-0970-9361
  affiliation: 3
- name: Téo Lemane
  orcid: 0000-0002-7210-3178
  affiliation: 1
- name: David Vallenet
  orcid: 0000-0001-6648-0332
  affiliation: 1
- name: Alexandra Calteau
  orcid: 0000-0002-5871-9347
  affiliation: 1
  corresponding: true

affiliations:
- name: LABGeM, Genomique Métabolique, CEA, Genoscope, Institut François Jacob, Université d’Évry, Université Paris-Saclay, CNRS, France.
  index: 1
  ror: 00xc55v17
- name: Laboratory of Biology and Modeling of the Cell, Ecole Normale Supérieure de Lyon, CNRS, Université Claude Bernard Lyon, Université de Lyon, France
  index: 2
  ror:  01bj4fd12
- name: Université Paris-Saclay, INRAE, MaIAGE, Jouy-en-Josas, France
  index: 3
  ror: 05qdnns64

date: 23 June 2026
bibliography: references.bib

---

# Summary

The exponential growth of genomic data, particularly for microbes, has made pangenomic approaches a gold standard for large-scale comparative genomics. By capturing the full genomic diversity of a species rather than relying on a single reference, pangenomics has transformed microbial genomics, revealing the adaptive potential of bacteria and the evolutionary dynamics underlying functional diversity. Among available tools, PPanGGOLiN distinguishes itself through its graph-based model coupled with statistical gene partitioning.

Here we present PPanGGOLiN v2, which introduces substantial improvements across three dimensions: new analytical features that expand what users can investigate, a comprehensive software architecture redesign that improves maintainability and extensibility, and performance improvements that address the computational demands of ever-growing genomic datasets.

# Statement of need

In the last decade, the number of sequenced microbial genomes has exploded from a few thousand to several millions. While this wealth of data offers immense potential, traditional genome-centric approaches reach their limits when exploring and interpreting variation at this scale [@land_insights_2015;@parks_standardized_2018]. The pangenome addresses this by condensing all genomic information from a set of genomes into a single structure, capturing not only gene presence/absence but also genetic context, e.g. mobile genetic elements or functional annotations, at the species or clade level [@tettelin_genome_2005].

PPanGGOLiN [@gautreau_ppanggolin_2020] addresses these needs through a graph-based model that clusters genes into families based on sequence similarity and encodes neighborhood information, enabling high compression of genomic diversity into a single data structure. Gene families are then partitioned using a Bernoulli Mixture Model coupled with a Markov Random Field that takes into account both family occurrence and genomic neighborhood. The software suite further integrates PanRGP [@bazin_panrgp_2020], for identification of Regions of Genomic Plasticity (RGPs) as well as their insertion spots, and PanModule [@bazin_panmodule_2021], for their description as conserved modules.

PPanGGOLiN v2 builds on this foundation by adding new analytical features: genomic context search, RGP clustering, pangenome projection and metadata association (\autoref{fig:overview}) together with a full architectural redesign and a modernized development workflow, making the tool more capable, maintainable, and open to community contributions.


# State of the field


Pangenomics has transformed prokaryotic genome analyses by moving beyond a mainly reference-centric approach toward new representations that capture the entire genomic diversity of microbial populations. Two representational models for pangenome graphs are widely used. Sequence-level graphs, with tools such as Bifrost [@holley_bifrost_2020] or the vg toolkit [@garrison_variation_2018], encode variation at the nucleotide level, providing resolution of any small variants, but at a high computational cost and interpretive complexity. Gene-level graphs with tools such as Panaroo [@tonkin-hill_producing_2020] or pangene [@li_exploring_2024] instead cluster genes into families to form graph nodes, with edges encoding genomic contiguity. This representation neglects fine-scale sequence variation but remains highly informative about gene presence/absence and mobile genetic elements at a fraction of the resource cost.

Most gene-level tools rely on a predefined occurrence threshold to split gene families into different partitions, which leads to inconsistent predictions depending on the quality and diversity of the available genomes. PPanGGOLiN addresses this using a statistical method relying on both the patterns of occurrence of gene families, and the pangenome graph topology to partition gene families into: *persistent*, *shell*, and *cloud* genomes. Among dozens of alternative gene-level pangenomic tools, only mOTUpan [@buck_motupan_2022] and micropan [@snipen_micropan_2015] use statistical methods to compute pangenome partitions.

Additionally, PPanGGOLiN builds on its graph-based representation and fine-grained partitioning to enable the identification of regions of genomic plasticity and of conserved modules — analyses that remain unique to PPanGGOLiN.


# New features

## Projection: pangenome annotation of external genomes

The projection feature enables annotation of external genomes using a previously computed pangenome, without recomputing the full pangenome structure. Genes from the input genome are assigned to existing gene families and partitions based on sequence similarity; genes without matches are considered genome-specific and assigned to the cloud partition. From these assignments, PPanGGOLiN further predicts Regions of Genomic Plasticity (RGPs), insertion spots, and conserved modules present in the projected genome. This is particularly useful for comparing newly sequenced genomes against an established reference pangenome.

In Version 2, PPanGGOLiN introduces a new output format compatible with Proksee [@grant_proksee_2023], enabling visualization of the circular genomes annotated together with pangenome information either from the pangenome itself or via projection allowing further interactive analyses directly within the Proksee platform.

## Genomic context extraction

A genomic context refers here to a set of genes conserved within similar genomic regions across multiple genomes. The pangenome graph already encodes neighborhood relationships and is therefore well suited for extracting such context. As input, PPanGGOLiN takes sequences corresponding to genes of interest, which are aligned to the reference sequence of each gene family to identify the target families in the pangenome. Once identified, gene families present within a defined genomic window around the target are extracted to construct a modified pangenome graph using a transitive closure taking into account a larger neighborhood.

The complete method is described in detail by Arnoux *et al* in the PANORAMA paper [@arnoux_panorama_2025]. This feature provides a species-level view of conserved genomic neighborhoods, enabling users to identify functionally related genes and, for example, accurately reconstruct metabolic pathways.


## RGP clustering

A new feature allows the comparative analysis of RGPs across genomes within the pangenome. Clustering is based on the similarity of gene content, where two RGPs are considered related if their genes belong to the same gene families. Similarity is quantified using a Gene Repertoire Relatedness (GRR) score.

Clustering is performed through graph modeling: each RGP is represented as a node, and an edge is added between two nodes if their GRR (min or max) exceeds a defined threshold (0.8 by default). Each connected component of the resulting graph corresponds to a cluster of RGPs sharing a common gene repertoire. Results are exportable as a tabulated file or in formats compatible with graph visualization tools such as Gephi[@bastian_gephi_2009] and Cytoscape [@shannon_cytoscape_2003]. This analysis facilitates the exploration of the diversity and dissemination of mobile genetic elements across genomes.

## Metadata

Another key addition is the ability to associate metadata with any pangenome element, including genes, contigs, genomes, gene families, edges, RGPs, spots, and modules. Metadata is provided as a tabulated file requiring at minimum the identifier of the target element, with no restriction on the type or number of additional fields, making the format highly flexible. Each metadata entry is linked to a named source, allowing the same element to carry multiple annotations from different sources simultaneously.

For performance and portability, metadata is stored directly in the pangenome file. While not used directly in pangenome computation, they are propagated to all PPanGGOLiN outputs, facilitating downstream interpretation and exploration of results.

Combined with RGP clustering, and pangenome projection, metadata association bridges the gap between pangenome structure and biological interpretation, allowing users to directly connect genomic variation to functional, ecological, or experimental data.

# Software design

PPanGGOLiN is built around a central data model in which the HDF5 pangenome file serves as the reference data object for the entire system. Analyses can be performed either step by step or through workflow commands that execute the pipeline in a single run, providing both flexibility and convenience (\autoref{fig:overview}). Once executed, the HDF5 pangenome file supports multiple uses: exporting to other file formats for downstream analyses, serving as input for additional PPanGGOLiN commands, and being accessed programmatically through the Python API in custom scripts. This central-object design also enhances reproducibility: parameters are resolved consistently (command line, then configuration file, then defaults), validated across analysis steps, and stored within the pangenome file, enabling analyses to be reliably reproduced and compared.

![**Overview of the PPanGGOLiN v2 workflow.** Each rounded box represents a command of the software. Commands marked with * are new in PPanGGOLiN v2. Commands are grouped into two sections. (*i*) *Pangenome construction*: starting from genomic data (with annotations or not), the `annotate`, `cluster`, `graph`, and `partition` commands progressively build the pangenome graph and partition gene families into *persistent*, *shell*, and *cloud* components. Then, two independent analyses can be performed: `rgp` followed by `spot` identify Regions of Genomic Plasticity (RGPs) and their insertion spots, while `module` identifies conserved genomic modules. All results are stored in a single HDF5 file (`pangenome.h5`), which serves as the central data object. Metadata from external sources can be associated with any pangenome element via the `metadata` command (\*), using a tabulated input file. (*ii*) *Pangenome exploitation*: the HDF5 file is used as input for four categories of downstream commands. *Pangenome analysis*: `msa` computes multiple sequence alignments for gene families; `rarefaction` computes the rarefaction curve of the pangenome; `rgp_cluster` (\*) clusters RGPs based on shared gene content. *Pangenome-guided analysis*: `align` maps an external genome or protein set to pangenome families; `context` (\*) extracts conserved genomic neighborhoods around genes of interest; `projection` (\*) annotates a new genome using an existing pangenome. *Pangenome output*: `draw`, `fasta`, `write_pangenome`, and `write_genomes` export results in various formats for downstream analyses and visualization.\label{fig:overview}](./ppanggolin_v2_overview.pdf)




# Technical enhancements

PPanGGOLiN was originally developed as a research project. Despite having robust and user-friendly workflows, its underlying codebase required a complete redesign and refactoring to meet long-term software engineering standards.

Several targeted improvements address performance and reproducibility. Prodigal [@hyatt_prodigal_2010] has been replaced by Pyrodigal [@larralde_pyrodigal_2022], reducing I/O overhead during genome annotation. The HDF5 pangenome file structure has been redesigned to drastically reduce file size (\autoref{fig:hdf5_size}) and improve programmatic access via API. Pangenome file reading was optimized, reducing loading time, especially for spots from over several minutes to few seconds for large datasets. Memory consumption during sequence writing was reduced by reading sequences directly from HDF5 files rather than loading them into memory. A configuration file system has been introduced to facilitate the integration of PPanGGOLiN into computational platforms (see Research impact), and error management has been overhauled to provide clearer and more informative user feedback.

Finally, documentation has been fully revised and is now hosted on ReadTheDocs, lowering the barrier for new users. PPanGGOLiN v2 is distributed via Bioconda[@gruning_bioconda_2018] and PyPI, ensuring broad accessibility across computing environments.

![**Pangenome HDF5 file size (GB) generated by PPanGGOLiN v1.2.74 vs v2.3.0** **across increasing genome counts.** PPanGGOLiN v2 consistently produces smaller files, supporting the impact of the v2 file structure redesign.\label{fig:hdf5_size}](./pangenome_h5_size.pdf)

# Research impact statement

PPanGGOLiN is widely adopted in microbial genomics, with \~300 citations on Google Scholar. Its results are integrated into MicroScope [@vallenet_microscope_2020], a microbial genome annotation and analysis platform, and a dedicated Galaxy [@the_galaxy_community_galaxy_2024] module has been developed to broaden accessibility. The tool serves as a core dependency for PANORAMA[@arnoux_panorama_2025], enabling cross-species biological system predictions and pangenome comparisons, and supports [PanGBank](https://pangbank.genoscope.cns.fr/) [@mainguy_pangbank_2026], a pangenome database. A companion project further extends this ecosystem by integrating PPanGGOLiN pangenomes — including RGPs, spots, and modules — alongside functional metadata annotations into a Neo4j graph database, enabling advanced graph-based queries across pangenome elements [@arnoux_integrating_2024]. PPanGGOLiN is actively maintained and evolves through community contributions, bug reports, and feature suggestions. In 2023, it received the Open Science Award for Open Source Research Software from the French Ministry of Higher Education and Research[^1].

[^1]: [https://www.ouvrirlascience.fr/the-2023-open-science-free-research-software-awards/](https://www.ouvrirlascience.fr/the-2023-open-science-free-research-software-awards/)

# AI usage disclosure

AI tools were used in a limited capacity during the later stages of development. In particular, GitHub Copilot (mainly using Claude Sonnet 4.6) was used for code review suggestions, assistance with writing some tests, and documentation. Claude Sonnet 4.6 (Anthropic) was used to assist in the design and iterative refinement of the workflow overview figure (\autoref{fig:overview}) and manuscript proofreading. The core software design, implementation, and architectural decisions were developed independently without AI assistance. All AI-generated outputs were systematically reviewed and validated by the authors.

# Acknowledgments

This research was supported in part by the CFR PhD program of the French Alternative Energies and Atomic Energy Commission (CEA) for JA and the BlueRemediomics project for JM, which is funded by the European Union under the Horizon Europe Program, Grant Agreement No. 101082304\.

We are grateful to the PPanGGOLiN user community for their contributions, feedback, and bug reports, which have continuously improved the tool and shaped the development of this new version.

# References
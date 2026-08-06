---
title: 'PPanGGOLiN V2: technical enhancement and extended functionalities for prokaryotic pangenome analysis'

tags:
- Python
- Bioinformatics
- Pangenomics
- Microbial genomics
- Mobile genetic elements
- Comparative genomics
- Pangenome graph

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
bibliography: paper.bib
---

# Summary

The exponential growth of genomic data, particularly for microbes, has made pangenomic approaches a gold standard for large-scale comparative genomics. By capturing the full genomic diversity of a species rather than relying on a single reference, pangenomics has transformed microbial genomics, revealing the adaptive potential of bacteria and the evolutionary dynamics underlying functional diversity. Among available tools, PPanGGOLiN distinguishes itself through its graph-based model coupled with statistical gene partitioning.

Here we present PPanGGOLiN v2, which introduces substantial improvements across three dimensions: new analytical features that expand what users can investigate, a comprehensive software architecture redesign that improves maintainability and extensibility, and performance improvements that address the computational demands of ever-growing genomic datasets.

# Statement of need

In the last decade, the number of sequenced microbial genomes has exploded from a few thousand to several millions. While this wealth of data offers immense potential, traditional genome-centric approaches reach their limits when exploring and interpreting variation at this scale [@land_insights_2015;@parks_standardized_2018]. A pangenome condenses all genomic information from a set of genomes into a single structure, capturing gene presence/absence and genetic context: e.g. mobile genetic elements or functional annotations, at the species or clade level [@tettelin_genome_2005].

PPanGGOLiN [@gautreau_ppanggolin_2020] addresses these needs through a graph-based model that clusters genes into families based on sequence similarity and encodes neighborhood information, compressing genomic diversity into a single structure. Gene families are then partitioned using a Bernoulli Mixture Model coupled with a Markov Random Field accounting for both family occurrence and genomic neighborhood. The software suite further integrates PanRGP [@bazin_panrgp_2020], for identification of Regions of Genomic Plasticity (RGPs) and their insertion spots, and PanModule [@bazin_panmodule_2021], for finding conserved modules.

PPanGGOLiN v2 builds on this foundation by adding new analytical features: genomic context search, RGP clustering, pangenome projection and metadata association (\autoref{fig:overview}) together with a full architectural redesign and a modernized development workflow, making the tool more capable, maintainable, and open to community contributions.

# State of the field

Pangenomics has transformed prokaryotic genome analysis by moving beyond a reference-centric approach toward representations that capture the entire genomic diversity of microbial populations. Two pangenome graphs models are widely used. Sequence-level graphs, like  Bifrost [@holley_bifrost_2020] or vg toolkit [@garrison_variation_2018], encode nucleotide-level variation, providing small variants resolution, but at high computational and interpretive complexity. Gene-level graphs, such as Panaroo [@tonkin-hill_producing_2020] or pangene [@li_exploring_2024], instead cluster genes into families to form graph nodes, with edges encoding genomic contiguity. This representation neglects fine-scale variation but remains highly informative about gene presence/absence and mobile elements at a fraction of the cost.

Most gene-level tools rely on a predefined occurrence threshold to split gene families into different partitions,leading to inconsistent predictions depending on  genome quality and diversity . PPanGGOLiN instead uses  a statistical method relying on both gene family occurrence patterns  and graph topology to partition families into: *persistent*, *shell*, and *cloud* genomes. Among dozens of gene-level pangenomic tools, only mOTUpan [@buck_motupan_2022] and micropan [@snipen_micropan_2015] also use statistical partitioning.

Additionally, PPanGGOLiN builds on its graph-based representation and fine-grained partitioning to identify RGPs, spots and  conserved modules; analyses unique to PPanGGOLiN.

# New features
## Projection: pangenome annotation of external genomes

The projection feature annotates external genomes using a previously computed pangenome, without recomputing the full structure. Input genes are assigned to existing gene families and partitions based on sequence similarity; unmatched genes are considered genome-specific and assigned to the cloud partition. From these assignments, PPanGGOLiN further predicts RGPs, spots, and modules present in the projected genome. This is useful for comparing newly sequenced genomes against an established reference pangenome.

In Version 2, PPanGGOLiN introduces a Proksee-compatible output format[@grant_proksee_2023],  to visualize the circular genomes annotated with pangenome information, from the pangenome itself or via projection, allowing further interactive analyses  within Proksee.

## Genomic context extraction

A genomic context refers to a set of genes conserved within similar genomic regions across multiple genomes. The pangenome graph already encodes neighborhood relationships, making it well suited for extracting such context. The complete method is described in detail by @arnoux_panorama_2026. This feature provides a species-level view of conserved genomic neighborhoods, letting users identify functionally related genes and, for example, accurately reconstruct metabolic pathways.

## RGP clustering

A new feature allows comparative analysis of RGPs across genomes in the pangenome. Clustering is based on gene content similarity, where two RGPs are considered related if their genes belong to the same families, quantified using a Gene Repertoire Relatedness (GRR) score.

Clustering is performed through graph modeling: each RGP is  a node, and an edge is added between two nodes if their GRR (min or max) exceeds a threshold (0.8 by default). Each connected component of the resulting graph corresponds to a cluster of RGPs sharing a common gene repertoire. Results are exportable as a tabulated file or in formats compatible with graph visualization tools such as Gephi[@bastian_gephi_2009] and Cytoscape [@shannon_cytoscape_2003], facilitating the exploration of the diversity and dissemination of mobile genetic elements across genomes.

## Metadata

Another key addition is the ability to associate metadata with any pangenome element, including genes, contigs, genomes, gene families, edges, RGPs, spots, and modules. Metadata is provided as a tabulated file requiring at minimum the target element’s identifier, with no restriction on the type or number of additional fields, making the format highly flexible. Each entry is linked to a named source, allowing an element to carry multiple annotations from different sources.

For performance and portability, metadata is stored directly in the pangenome file. While not used directly in pangenome computation, it is propagated to all outputs, facilitating downstream interpretation and exploration.

Combined with RGP clustering, and pangenome projection, metadata association bridges the gap between pangenome structure and biological interpretation, letting users connect genomic variation to functional, ecological, or experimental data.

# Software design

PPanGGOLiN is built around a central data model in which the HDF5 pangenome file serves as the reference data object for the entire system. Analyses can be performed step by step or through workflow commands executing the pipeline in a single run, providing flexibility and convenience (\autoref{fig:overview}). Once executed, the HDF5 pangenome file supports multiple uses: export to other formats for downstream analyses, input for additional PPanGGOLiN commands, and programmatic access through the Python API. This central-object design also enhances reproducibility: parameters are resolved consistently (command line, then configuration file, then defaults), validated across analysis steps, and stored within the pangenome file, enabling analyses to be reliably reproduced and compared.

![**Overview of the PPanGGOLiN v2 workflow.** Rounded boxes represent software commands. Commands marked with * are new in v2. Commands fall into two sections. (*i*) *Pangenome construction*: from genomic data (annotated or not),  to progressively build the pangenome graph and partition gene families into *persistent*, *shell*, and *cloud*. Two independent analyses then follow: `rgp` and `spot` identify RGPs and their insertion spots, while `module` identifies conserved genomic modules. Results are stored in a single HDF5 file (`pangenome.h5`), the central data object. External metadata can be associated with any pangenome element via  `metadata` (\*), using a tabulated file. (*ii*) *Pangenome exploitation*: the HDF5 file serves as input for four categories of downstream commands. *Pangenome analysis*: `msa` computes multiple sequence alignments for gene families; `rarefaction` computes the pangenome’s rarefaction curve; `rgp_cluster` (\*) clusters RGP based on shared gene content. *Pangenome-guided analysis*: `align` maps an external genome or protein set to pangenome families; `context` (\*) extracts conserved genomic neighborhoods around genes of interest; `projection` (\*) annotates a new genome using an existing pangenome. *Pangenome output*: `draw`, `fasta`, `write_pangenome`, and `write_genomes` export results in various formats for downstream analyses and visualization.\label{fig:overview}](./ppanggolin_v2_overview.pdf)

# Technical enhancements

PPanGGOLiN was originally developed as a research project. Despite robust and user-friendly workflows, its codebase required a complete redesign to meet long-term software engineering standards.

Several targeted improvements address performance and reproducibility. Prodigal [@hyatt_prodigal_2010] was replaced by Pyrodigal [@larralde_pyrodigal_2022], reducing I/O overhead during annotation. The HDF5 pangenome file structure was redesigned to drastically reduce file size (\autoref{fig:hdf5_size}) and improve API access. Pangenome file reading was optimized, reducing loading time from minutes to seconds for large datasets. Memory consumption during sequence writing was reduced by reading sequences directly from HDF5 files rather than loading them into memory. A configuration file system now facilitates PPanGGOLiN  integration into computational platforms (see Research impact).

Finally, documentation was fully revised and is now hosted on ReadTheDocs, lowering the barrier for new users. PPanGGOLiN v2 is distributed via Bioconda[@gruning_bioconda_2018] and PyPI, ensuring broad accessibility across computing environments.

![**Pangenome HDF5 file size (GB) generated by PPanGGOLiN v1.2.74 vs v2.3.0** **across increasing genome counts.** PPanGGOLiN v2 consistently produces smaller files, supporting the impact of the v2 file structure redesign.\label{fig:hdf5_size}](./pangenome_h5_size.pdf)

# Research impact statement

PPanGGOLiN is widely adopted in microbial genomics, with \~300 citations on Google Scholar. Its results are integrated into MicroScope [@vallenet_microscope_2020], a microbial genome annotation and analysis platform, and a dedicated Galaxy [@the_galaxy_community_galaxy_2024] module  broadens accessibility. The tool  is a core dependency for PANORAMA [@arnoux_panorama_2026], enabling cross-species biological system predictions and pangenome comparisons, and supports [PanGBank](https://pangbank.genoscope.cns.fr/), a pangenome database. A companion project further extends this ecosystem by integrating PPanGGOLiN pangenomes — including RGP, spots, and modules — alongside functional metadata into a Neo4j graph database, enabling advanced graph-based queries across pangenome elements [@arnoux_integrating_2024]. PPanGGOLiN is actively maintained and evolves through community contributions, bug reports, and feature suggestions. In 2023, it received the [Open Science Award for Open Source Research Software](https://www.ouvrirlascience.fr/the-2023-open-science-free-research-software-awards/) from the French Ministry of Higher Education and Research.

# AI usage disclosure

AI tools were used in a limited capacity during later development. GitHub Copilot (mainly Claude Sonnet 4.6) assisted with code review suggestions, test writing, and documentation. Claude Sonnet 4.6 (Anthropic) assisted in designing and refining the workflow overview figure (\autoref{fig:overview}) and in manuscript proofreading. Core software design, implementation, and architectural decisions were developed independently without AI assistance. All AI-generated outputs were systematically reviewed and validated by authors.

# Acknowledgments

This research was supported by the CFR PhD program of the French Alternative Energies and Atomic Energy Commission (CEA) for JA and the BlueRemediomics project for JM, funded by the European Union under the Horizon Europe Program, Grant Agreement No. 101082304\.

We thank the PPanGGOLiN user community for their contributions, feedback, and bug reports, which have continuously improved the tool and shaped this new version.

# References

# Projection

The ppanggolin projection command allows you to annotate external genomes using an existing pangenome. This process eliminates the need to recompute all components, streamlining the annotation process. Input genomes are expected to belong to the same species.

Genes within the input genome are aligned with genes in the pangenome to determine their gene families and partitions. Genes that do not align with any existing gene in the pangenome are considered specific to the input genome and are assigned to the "Cloud" partition. Based on the alignment and partition assignment, Regions of Plasticity (RGPs) within the input genome are predicted. Each RGP that is not located on a contig border is assigned to a spot of insertion. Finally, conserved modules of the pangenome found in the input genome are reported in the output files.

## Input files:

This command supports two input modes depending on whether you want to project a single genome or multiple genomes at once:

Multiple Files in One TSV:
- **Options**: `--fasta` or `--anno`
- **Description**: You can provide a tab-separated file listing organism names alongside their respective FASTA genomic sequences or annotation filepaths, with one line per organism. This mode is suitable when you want to annotate multiple genomes in a single operation. The format of this file is identical to the format used in the annotate and workflow commands; for more details, refer here.

Single File:
- **Options**: `--organism_name` with `--fasta` or `--anno` and `--circular_contigs` (optional)
- **Description**: When annotating a single genome, you can directly provide a single FASTA genomic sequence file or an annotation file in GFF/GBFF format. Additionally, specify the name of the organism using the `--organism_name` option. You can also indicate circular contigs using the `--circular_contigs` option when necessary.


## Output files:

The Output directory contains `summary_projection.tsv` giving an overview of the projection. one line per organism.


| Column                               | Description|
|--------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Organism name                        | This column contains name or identifier of the organisms being analyzed.|
| Pangenome file                       | The path to the pangenome file (pangenome.h5) used for the analysis.|
| Contigs                              | The number of contigs in the projected genome.|
| Genes                                | The total number of genes identified in the input genome.|
| Families                             | The total number of gene families to which genes in the genome of the input organism are assigned.|
| Persistent genes                     | The number of genes in the "Persistent" partition.|
| Persistent families                  | The number of gene families in the "Persistent" partition.|
| Shell genes                          | The number of genes in the "Shell" partition.|
| Shell families                       | The number of gene families in the "Shell" partition.|
| Cloud genes                          | The number of genes in the "Cloud" partition.|
| Cloud families                       | The number of gene families in the "Cloud" parition.|
| Cloud specific families              | The number of gene families that are specific to the input organism. These families are unique to the input organism and do not have homologs in any other genomes within the pangenome and have been assigned to the "Cloud" partition.|
| completeness           | This indicates the proportion of single copy markers from the persistent partition that are present in the genome. While it is expected to be relatively close to 100 when working with isolates, it may be particularly interesting when working with very fragmented genomes as this provide a *de novo* estimation of the completeness based on the expectation that single copy markers within the persistent should be mostly present in all individuals of the studied taxonomic group. |
| RGPs (Regions of Genomic Plasticity) | The number of Regions of Genomic Plasticity (RGPs) predicted within the input genome.|
| Spots                                | The total number of spots of insertion associated with RGPs in the input genome.|
| New spots                            | The number of new insertion spots that have been identified in the input genome. These spots represent novel genomic regions compared to other genomes in the pangenome.|
| Modules                              | The number of modules that have been projected onto the input genome.|


Additionally, within the Output directory, there is a subdirectory for each input genome, named after the input genome itself. Each of these subdirectories contains several files:

For Gene Family and Partition of Input Genes:

- `cds_sequences.fasta`: This file contains the sequences of coding regions (CDS) from the input genome.
- `gene_to_gene_family.tsv`: It provides the mapping of genes to gene families of the pangenome. its format follows [this output](Outputs.md#gene-families-and-genes)
- `sequences_partition_projection.tsv`: This file maps the input genes to its partition (Persistent, Shell or Cloud).
- `specific_genes.tsv`: This file list the gene of the input genomes that do not align to any gene of the pangenome. These genes are assigned to Cloud parititon. 

For RGPs and Spots:

- `plastic_regions.tsv`: This file contains information about Regions of Genomic Plasticity (RGPs) within the input genome. Its format follows [this output](Outputs.md#plastic-regions).
- `input_organism_rgp_to_spot.tsv`: It provides information about the association between RGPs and insertion spots in the input genome. Its format follows [this ouput](Outputs.md#spots).

Optionally, you can produce a graph of the RGPs using the `--spot_graph` option. This graph is similar as the one produce by the `ppanggolin spot` command.

For Modules:

- `modules_in_input_organism.tsv`: This file lists the modules that have been found in the input genome. Its format follows [this ouput](Outputs.md#modules-in-organisms).


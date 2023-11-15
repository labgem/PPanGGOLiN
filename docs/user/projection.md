# Projection
The ppanggolin projection command allows you to annotate external genomes using an existing pangenome. This process eliminates the need to recompute all components, streamlining the annotation process. Input genomes are expected to belong to the same species.

Genes within the input genome are aligned with genes in the pangenome to determine their gene families and partitions. Genes that do not align with any existing gene in the pangenome are considered specific to the input genome and are assigned to the "Cloud" partition. The number of this specific cloud families are detailed in the summary table. The summary table provides a count of these specific cloud families. 

Based on the alignment and partition assignment, Regions of Plasticity (RGPs) within the input genome are predicted. Each RGP that is not located on a contig border is assigned to a spot of insertion. Finally, conserved modules of the pangenome found in the input genome are reported in the output files.

## Input files:

This command supports two input modes depending on whether you want to project a single genome or multiple genomes at once:

Multiple Files in One TSV:
- **Options**: `--fasta` or `--anno`
- **Description**: You can provide a tab-separated file listing organism names alongside their respective FASTA genomic sequences or annotation filepaths, with one line per organism. This mode is suitable when you want to annotate multiple genomes in a single operation. The format of this file is identical to the format used in the annotate and workflow commands; for more details, refer here.

Single File:
- **Options**: `--organism_name` with `--fasta` or `--anno` and `--circular_contigs` (optional)
- **Description**: When annotating a single genome, you can directly provide a single FASTA genomic sequence file or an annotation file in GFF/GBFF format. Additionally, specify the name of the organism using the `--organism_name` option. You can also indicate circular contigs using the `--circular_contigs` option when necessary.


## Output Files

Within the Output directory, the `summary_projection.tsv` file provides an overview of the projection, featuring one line per organism. This file includes all the columns described in the [organisms-statistics](Outputs.md#organisms-statistics) output section, along with specific projection-related columns detailed below:

| Column                      | Description                                                                                   |
|-----------------------------|-----------------------------------------------------------------------------------------------|
| Cloud specific families     | Number of cloud-specific families. These gene families do not match any existing families within the pangenome. |
| Pangenome file              | The file path of the pangenome used in the projection.                                           |
| New spots                   | Number of newly identified spots in the input genome.                                            |



The `summary_projection.tsv` file contains additional columns pertinent to the projection process, supplementing the standard organism statistics with specialized information about the projection, facilitating a comprehensive understanding of the analyzed data.

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


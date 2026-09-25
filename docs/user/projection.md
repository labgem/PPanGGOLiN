# Projection

The ppanggolin projection command allows you to annotate external genomes using an existing pangenome. This process eliminates the need to recompute all components, streamlining the annotation process. Input genomes are expected to belong to the same species.

Genes within the input genome are aligned with genes in the pangenome to determine their gene families and partitions. Genes that do not align with any existing gene in the pangenome are considered specific to the input genome and are assigned to the "Cloud" partition. The numbers of these specific cloud families are detailed in the summary table.

Based on the alignment and partition assignment, Regions of Plasticity (RGPs) within the input genome are predicted. Each RGP that is not located on a contig border is assigned to a spot of insertion. Finally, conserved modules of the pangenome found in the input genome are reported in the output files.


### Input files

Projection accepts the same genome input formats as the `annotate` and `all` commands. Genomes can be provided either as genomic FASTA files using `--fasta` or as GFF3/GBFF annotation files using `--anno`.

You can project either a single genome or multiple genomes at once.

### Project multiple genomes

To project multiple genomes, provide `--fasta` or `--anno` with a tab-separated file listing one genome per line. The file format is identical to the input files used to create a pangenome.

Examples with 50 *Chlamydia trachomatis* genomes are available in the [testingDataset](https://github.com/labgem/PPanGGOLiN/tree/master/testingDataset) directory: [`genomes.gbff.list`](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/genomes.gbff.list) for annotation files and [`genomes.fasta.list`](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/genomes.fasta.list) for FASTA files.

```bash
ppanggolin projection -p pangenome.h5 --anno genomes.gbff.list
```

or:

```bash
ppanggolin projection -p pangenome.h5 --fasta genomes.fasta.list
```

### Project a single genome

For a single genome, provide the file directly and specify its name with `--genome_name`. The input can be either a genomic FASTA file or a GFF3/GBFF annotation file.

```bash
ppanggolin projection -p pangenome.h5 \
    --fasta genome_A.fasta \
    --genome_name genome_A
```

or:

```bash
ppanggolin projection -p pangenome.h5 \
    --anno genome_A.gbff \
    --genome_name genome_A
```

The optional `--circular_contigs` argument can be used to specify circular contigs.

## Output Files

Within the Output directory, the `summary_projection.tsv` file provides an overview of the projection, featuring one line per genome. This file includes all the columns described in the [genome-statistics table](./PangenomeAnalyses/pangenomeStat.md#genome-statistics-table) section, along with specific projection-related columns detailed below:

| Column                  | Description                                                                                                     |
|-------------------------|-----------------------------------------------------------------------------------------------------------------|
| Cloud_specific_families | Number of cloud-specific families. These gene families do not match any existing families within the pangenome. |
| Pangenome_file          | The file path of the pangenome used in the projection.                                                          |
| New_spots               | Number of newly identified spots in the input genome.                                                           |



Additionally, within the Output directory, there is a subdirectory for each input genome, named after the input genome itself. Each of these subdirectories contains several files:


For Gene Family and Partition of Input Genes:

- `cds_sequences.fasta`: This file contains the sequences of coding regions (CDS) from the input genome.
- `gene_to_gene_family.tsv`: It provides the mapping of genes to gene families of the pangenome. Its format follows [this output](PangenomeAnalyses/pangenomeStat.md#gene-families-to-gene-associations).
- `sequences_partition_projection.tsv`: This file maps the input genes to its partition (Persistent, Shell or Cloud).
- `specific_genes.tsv`: This file lists the gene of the input genomes that do not align to any gene of the pangenome. These genes are assigned to Cloud partition. 

For RGPs and Spots:

- `plastic_regions.tsv`: This file contains information about RGPs within the input genome. Its format follows [this output](RGP/rgpOutputs.md#rgp-outputs).
- `input_genome_rgp_to_spot.tsv`: It provides information about the association between RGPs and insertion spots in the input genome. Its format follows [this output](RGP/rgpOutputs.md#summarize-spots).

Optionally, you can generate a graph of the spots using the `--spot_graph` option. This graph resembles the one produced by the `ppanggolin draw --spots` command, which is detailed [here](RGP/rgpOutputs.md#draw-spots).

For Modules:

- `modules_in_input_genome.tsv`: This file lists the modules that have been found in the input genome. Its format follows [this output](Modules/moduleOutputs.md#module-outputs).




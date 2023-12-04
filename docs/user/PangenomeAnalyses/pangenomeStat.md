### Pangenome statistics

The command `ppanggolin write_pangenome` allows to write 'flat' files that describe the pangenome and its elements.

#### Organism statistics table

The `organisms_statistics.tsv` file is a tab-separated file summarizing the content of each of the genomes used for building the pangenome. This file is useful when working with fragmented data, such as MAGs, or when investigating potential outliers within your dataset, such as chimeric or taxonomically disparate genomes.

The first lines starting with a `#` are indicators of parameters used when generating the numbers describing each organisms, and should not be read if loading this into a spreadsheet. They will be skipped automatically if you load this file with R.

This file comprises 32 columns described in the following table:

| Column                      | Description                                                                                   |
|-----------------------------|-----------------------------------------------------------------------------------------------|
| Organism_name               | Name of the organism to which the provided genome belongs                                     |
| Contigs                     | Number of contigs present in the genome                                                       |
| Genes                       | Total number of genes in the genome                                                           |
| Fragmented_genes            | Number of genes flagged as fragmented. Refer to the [defragmentation](./pangenomeCluster.md#defragmentation) section for detailed information on the fragmentation process.    |
| Families                    | Total number of gene families present in the genome                                           |
| Families_with_fragments     | Number of families containing fragmented genes                                                |
| Families_in_multicopy       | Number of families present in multiple copies                                                 |
| Soft_core_families          | Number of families categorized as soft core                                                   |
| Soft_core_genes             | Number of genes within soft core families                                                     |
| Exact_core_families         | Number of families categorized as exact core                                                  |
| Exact_core_genes            | Number of genes within exact core families                                                    |
| Persistent_genes            | Number of genes classified as persistent                                                      |
| Persistent_fragmented_genes | Number of genes flagged as persistent and fragmented                                          |
| Persistent_families         | Number of families categorized as persistent                                                  |
| Persistent_families_with_fragments | Number of persistent families containing fragmented genes                              |
| Persistent_families_in_multicopy | Number of persistent families present in multiple copies                                 |
| Shell_genes                 | Number of genes classified as shell                                                           |
| Shell_fragmented_genes      | Number of genes flagged as shell and fragmented                                               |
| Shell_families              | Number of families categorized as shell                                                       |
| Shell_families_with_fragments | Number of shell families containing fragmented genes                                        |
| Shell_families_in_multicopy | Number of shell families present in multiple copies                                           |
| Cloud_genes                 | Number of genes classified as cloud                                                           |
| Cloud_fragmented_genes      | Number of genes flagged as cloud and fragmented                                               |
| Cloud_families              | Number of families categorized as cloud                                                       |
| Cloud_families_with_fragments | Number of cloud families containing fragmented genes                                        |
| Cloud_families_in_multicopy | Number of cloud families present in multiple copies                                           |
| Completeness                | Proportion of persistent families present in the genome; expected to be close to 100 for isolates |
| Contamination               | Proportion of single copy persistent families found in multiple copy in the genome.  |
| Fragmentation               | Proportion of families with fragmented genes in the genome |
| RGPs                        | Number of Regions of Genomic Plasticity identified                                            |
| Spots                       | Number of spot IDs in which the RGPs are inserted                                             |
| Modules                     | Number of module IDs within gene families                                  |




```{note}
If you have predicted RGPs, spots or modules in your pangenome, corresponding columns will be added.
```


This table can be generated using the `write_pangenome` subcommand with the flag `--stats` as such : 

```bash
ppanggolin write_pangenome -p pangenome.h5 --stats
```


The flag `--stats` will also generate the `mean_persistent_duplication.tsv` file desdcribe [here](#mean-persistent-duplication).



##### Genome Metrics Overview

**- Completeness**
The completeness value is expected to be relatively close to 100 when working with isolates, it may be particularly interesting when working with very fragmented genomes as this provides a *de novo* estimation of the completeness based on the expectation that persistent genes should be mostly present in all individuals of the studied taxonomic group

**- Contamination**
The Contamination value represents the proportion of single-copy persistent families found in multiple copies within the genome. These single-copy persistent families are computed using the `duplication_margin` parameter specified at the beginning of the file. They encompass all persistent gene families present in single copy in less than 5% of the genomes by default. 

In this computation, fragmented genes are excluded. Therefore, if a family exists in multicopy due to fragmented genes, it will still be counted as single copy. Contamination assessment is particularly useful for identifying potential chimeric genomes, especially in lower-quality genomes.

**- Fragmentation**
The fragmentation value denotes the proportion of families containing fragmented genes within the genome. A high fragmentation value may indicate a highly fragmented genome.





#### Mean Persistent Duplication

The `mean_persistent_duplication.tsv` is a tab-separated file that lists the gene families and their duplication ratio, their mean presence in the pangenome and whether it is considered a 'single copy marker'. A gene family is considered duplicated it is found in single copy in less than 5% of the genomes by default. This threshold can be adjusted with the parameter `--dup_margin`. And the value that has been used to generated this file is specififed as a comment line strating with a '#'. This notion of single copy marker is used to compute a contamination value for each genome in the [organisms statistics table](#organism-statistics-table) described previously, where the contamination is the proportion of single copy marker found in multicopy in a specified genome. 




To generate this file, use the following `write_pangenome` subcommand:

```bash
ppanggolin write_pangenome -p pangenome.h5 --stats
```

Executing this command will also create the `organisms_statistics.tsv` file, detailed in the section labeled [here](#mean-persistent-duplication).


### Mean Persistent Duplication

The `mean_persistent_duplication.tsv` file lists the gene families along with their duplication ratios, average presence in the pangenome, and classification as 'single copy markers.' In this context, a gene family is not considered in single copy if it appears in single copy in less than 5% of the genomes by default. This default threshold can be adjusted using the `--dup_margin` parameter. The chosen threshold value for generating this file is indicated within a comment line starting with a '#'.

This notion of single copy markers is used for calculating contamination values in the [organisms statistics table](#organism-statistics-table) described earlier.

Below an example excerpt from this file:

```tsv
#duplication_margin=0.05
persistent_family	duplication_ratio	mean_presence	is_single_copy_marker
J4H57_RS02250	0.0	1.0	True
K6U54_RS13115	0.0	1.0	True
J4H33_RS19875	0.0	1.0	True
J4H71_RS19770	0.0	1.0	True
NB568_RS20780	0.0	1.0	True
JCM18904_4793	0.0	1.0	True
JCM18904_4607	0.0	1.0	True
JCM18904_2671	0.003	1.003	True
JCM18904_2527	0.0	1.0	True
K04M1_RS08450	0.698	1.698	False
[...]

```

The `mean_persistent_duplication.tsv` file, can be generated using the `write_pangenome` subcommand with the flag `--stats` as such : 

```bash
ppanggolin write_pangenome -p pangenome.h5 --stats
```


The flag `--stats` will also generate the `organisms_statistics.tsv` file desdcribe [here](#organism-statistics-table).



(gene-presence-absence)=
### Gene Presence-Absence Matrix

The `gene_presence_absence.Rtab` file represents a presence-absence matrix wherein columns are the genomes used to construct the pangenome, and rows correspond to gene families. Each gene family is identified by the identifier of their representative gene.

The matrix contains '1' if the gene family is present in a particular genome and '0' if absent. This format mirrors the structure of the `gene_presence_absence.Rtab` file obtained from the pangenome analysis software [Roary](https://sanger-pathogens.github.io/Roary/).

To generate this file, use the `write_pangenome` subcommand with the `--Rtab` flag as follows:

```bash
ppanggolin write_pangenome -p pangenome.h5 --Rtab
```


### Matrix File
The `matrix.csv` file, formatted as a .csv file, follows a structure similar to the `gene_presence_absence.csv` file generated by [Roary](https://sanger-pathogens.github.io/Roary/). This file format is compatible with [Scoary](https://github.com/AdmiralenOla/Scoary) for performing pangenome-wide association studies.

To generate this file, use the `write_pangenome` subcommand with the `--csv` flag:

```bash
ppanggolin write_pangenome -p pangenome.h5 --csv
```



### Partitions Files

The 'Partitions' files are stored within the `partitions` directory and are named after the specific partition they represent (e.g., 'persistent.txt' for the persistent partition). Each file contains a list of gene family identifiers corresponding to the gene families belonging to that particular partition. The format consists of one family identifier per line, facilitating their usage in upstream analysis workflows.

To generate these files, use the `write_pangenome` subcommand with the `--partitions` flag:

`ppanggolin write_pangenome -p pangenome.h5 --partitions`


### Gene Families to Genes Associations

The `gene_families.tsv` file mirrors the format provided by [MMseqs2](https://github.com/soedinglab/MMseqs2) through its `createtsv` subcommand. This file structure comprises three columns: the gene family name in the first column, the gene names in the second, and a third column that remains empty or contains an "F" to denote potential gene fragments instead of complete genes. This indication appears only if the [defragmentation](./pangenomeCluster.md#defragmentation) pipeline has been used.

To generate this file, use the `write_pangenome` subcommand with the `--families_tsv` flag:

`ppanggolin write_pangenome -p pangenome.h5 --families_tsv`

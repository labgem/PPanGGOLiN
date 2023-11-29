### Pangenome statistics


The command `ppanggolin write_pangenome` allows to write 'flat' files that describe the pangenome and its elements. 

#### Organism statistics table

The 'organisms_statistics.tsv' file is a tab-separated file summarizing the content of each of the genomes used for building the pangenome. This file is useful when working with fragmented data, such as MAGs, or when investigating potential outliers within your dataset, such as chimeric or taxonomically disparate genomes.

The first lines starting with a '#' are indicators of parameters used when generating the numbers describing each organisms, and should not be read if loading this into a spreadsheet. They will be skipped automatically if you load this file with R.

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
If you already have predict RGPs, spots or modules in your pangenome, a corresponding column will be added.
```


It can be generated using the 'write_pangenome' subcommand as such : 

```bash
ppanggolin write_pangenome -p pangenome.h5 --stats
```


This command will also generate the 'mean_persistent_duplication.tsv' file as desdcribe [here](#mean-persistent-duplication).



##### Genome Metrics Overview

**- Completeness**
The completeness value is expected to be relatively close to 100 when working with isolates, it may be particularly interesting when working with very fragmented genomes as this provides a *de novo* estimation of the completeness based on the expectation that persistent genes should be mostly present in all individuals of the studied taxonomic group

**- Contamination**
The Contamination value represents the proportion of single-copy persistent families found in multiple copies within the genome. These single-copy persistent families are computed using the `duplication_margin` parameter specified at the beginning of the file. They encompass all persistent gene families present in single copy in less than 5% of the genomes by default. 

In this computation, fragmented genes are excluded. Therefore, if a family exists in multicopy due to fragmented genes, it will still be counted as single copy. Contamination assessment is particularly useful for identifying potential chimeric genomes, especially in lower-quality genomes.

**- Fragmentation**
The fragmentation value denotes the proportion of families containing fragmented genes within the genome. A high fragmentation value may indicate a highly fragmented genome.




#### Mean persistent duplication

The file `mean_persistent_duplication.tsv` is a  tab-separated file, with a single parameter written as a comment at the beginning of the file, which indicates the proportion of genomes in which a gene family must be present more than once to be considered 'duplicated' (and not single copy marker).
This file lists the gene families, their duplication ratio, their mean presence in the pangenome and whether it is considered a 'single copy marker' or not, which is particularly useful when calculating the contamination recorded in the [organisms statistics table](#organism-statistics-table) described previously.

It can be generated using the 'write_pangenome' subcommand : 

```bash
ppanggolin write_pangenome -p pangenome.h5 --stats
```

This command will also generate the 'organisms_statistics.tsv' file as desdcribe [here](#mean-persistent-duplication).

(gene-presence-absence)=
#### Gene presence absence

This file is basically a presence absence matrix. The columns are the genomes used to build the pangenome, the lines are the gene families. The identifier of the gene family is the gene identifier chosen as a representative.
 There is a 1 if the gene family is present in a genome, and 0 otherwise. It follows the exact same format than the 'gene_presence_absence.Rtab' file that you get from the pangenome analysis software [Roary](https://sanger-pathogens.github.io/Roary/)

It can be generated using the 'write' subcommand as such: 

`ppanggolin write_pangenome -p pangenome.h5 --Rtab`

#### Matrix

This file is a .csv file following a format alike the  gene_presence_absence.csv file generated by [roary](https://sanger-pathogens.github.io/Roary/), and works with [scoary](https://github.com/AdmiralenOla/Scoary) if you want to do pangenome-wide association studies.

It can be generated using the 'write' subcommand as such: 

`ppanggolin write_pangenome -p pangenome.h5 --csv`


#### Partitions

Those files will be stored in the 'partitions' directory and will be named after the partition that they represent (like persistent.txt for the persistent partition). In each of those file there will be a list of gene family identifiers that correspond to the gene families belonging to that partition, one family per line, should you need it for your pipelines or during your analysis.

You can generate those files as such:  

` ppanggolin write_pangenome -p pangenome.h5 --partitions`


#### Gene families to genes associtations

You can write a list containing the gene family assigned to every single gene of your pangenome, in a file format extactly like the one provided by [MMseqs2](https://github.com/soedinglab/MMseqs2) through its subcommand 'createtsv'.
It is basically a three-column file listing the gene family name in the first column, and the gene names in the second. A third column is either empty, or has an "F" in it. In that case it indicates that the gene is potentially a gene fragment and not complete. This will be indicated only if the [defragmentation](./pangenomeCluster.md#defragmentation) pipeline is used.

You can obtain it as such:  

`ppanggolin write_pangenome -p pangenome.h5 --families_tsv`
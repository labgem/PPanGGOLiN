### Pangenome statistics and metrics
#### Organism Statistique

The organisms_statistics.tsv file is a tab-separated file describing the content of each of the genome used for building the pangenome. It might be useful when working with fragmented data such as MAGs or if you suspect some of your genomes to be chimeric, or to not belong to your taxonomic group (as those genomes will be outliers regarding to the numbers in this file).
The first lines starting with a '#' are indicators of parameters used when generating the numbers describing each organisms, and should not be read if loading this into a spreadsheet. They will be skipped automatically if you load this file with R.

This file is made of 15 columns described in the following table

| Column                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| organism               | Indicates the organism's name to whom the provided genome belongs to                                                                                                                                                                                                                                                                                                                                                                                                                    |
| nb_families            | Indicates the number of gene families present in that genome                                                                                                                                                                                                                                                                                                                                                                                                                            |
| nb_persistent_families | The number of persistent families present in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                |
| nb_shell_families      | The number of shell families present in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| nb_cloud_families      | The number of cloud families present in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| nb_exact_core          | The number of exact core families present in that genome. This number should be identical in all genomes.                                                                                                                                                                                                                                                                                                                                                                               |
| nb_soft_core           | The number of soft core families present in that genome. The threshold used is indicated in the #soft_core line at the beginning of the file, and is 0.95 by default.                                                                                                                                                                                                                                                                                                                   |
| nb_genes               | The number of genes in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| nb_persistent_genes    | The number of genes whose family is persistent in that genome                                                                                                                                                                                                                                                                                                                                                                                                                           |
| nb_shell_genes         | The number of genes whose family is shell in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                |
| nb_cloud_genes         | The number of genes whose family is cloud in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                |
| nb_exact_core_genes    | The number of genes whose family is exact core in that genome                                                                                                                                                                                                                                                                                                                                                                                                                           |
| nb_soft_core_genes     | The number of genes whose family is soft core in that genome                                                                                                                                                                                                                                                                                                                                                                                                                            |
| completeness           | This is an indicator of the proportion of single copy markers in the persistent that are present in the genome. While it is expected to be relatively close to 100 when working with isolates, it may be particularly interesting when working with very fragmented genomes as this provide a *de novo* estimation of the completess based on the expectation that single copy markers within the persistent should be mostly present in all individuals of the studied taxonomic group |
| nb_single_copy_markers | This indicates the number of present single copy markers in the genomes. They are computed using the parameter duplication_margin indicated at the beginning of the file. They correspond to all of the persistent gene families that are not present in more than one copy in 5% (or more) of the genomes by default.                                                                                                                                                                  |

It can be generated using the 'write' subcommand as such: 

`ppanggolin write -p pangenome.h5 --stats`

This command will also generate the 'mean_persistent_duplication.tsv' file as desdcribe [here](#mean-persistent-duplication).

```{note}
If you already have predict RGPs, spots or modules in your pangenome, a corresponding column will be added.
```

#### Mean persistent duplication
This file is a .tsv file, with a single parameter written as a comment at the beginning of the file, which indicates the proportion of genomes in which a gene family must be present more than once to be considered 'duplicated' (and not single copy marker).
This file lists the gene families, their duplication ratio, their mean presence in the pangenome and whether it is considered a 'single copy marker' or not, which is particularly useful when calculating the completeness recorded in the [organisms statistics](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#organisms-statistics) file described previously.


#### Genomic fluidity
The genomic fluidity is described as *a robust metric to categorize the
gene-level similarity among groups of sequenced isolates.* 
[more information here](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-32)

We add the possibility to get genomic fluidity for all the pangenome or 
for specific partition. The genomic fluidity is computable like that (with below an example result):
```
ppanggolin metrics -p pangenome --genome_fluidity
...
Genomes fluidity: all=0.026, shell=0.477, cloud=0.045, accessory=0.554
```
*all* correspond to all the family in the pangenome (core and accessory)

(gene-presence-absence)=
#### Gene presence absence

This file is basically a presence absence matrix. The columns are the genomes used to build the pangenome, the lines are the gene families. The identifier of the gene family is the gene identifier chosen as a representative.
 There is a 1 if the gene family is present in a genome, and 0 otherwise. It follows the exact same format than the 'gene_presence_absence.Rtab' file that you get from the pangenome analysis software [roary](https://sanger-pathogens.github.io/Roary/)

It can be generated using the 'write' subcommand as such: 

`ppanggolin write -p pangenome.h5 --Rtab`

## Matrix

This file is a .csv file following a format alike the  gene_presence_absence.csv file generated by [roary](https://sanger-pathogens.github.io/Roary/), and works with [scoary](https://github.com/AdmiralenOla/Scoary) if you want to do pangenome-wide association studies.

It can be generated using the 'write' subcommand as such: 

`ppanggolin write -p pangenome.h5 --csv`


#### Partitions

Those files will be stored in the 'partitions' directory and will be named after the partition that they represent (like persistent.txt for the persistent partition). In each of those file there will be a list of gene family identifiers that correspond to the gene families belonging to that partition, one family per line, should you need it for your pipelines or during your analysis.

You can generate those files as such:  

` ppanggolin write -p pangenome.h5 --partitions`

## Projection

This option writes in a 'projection' directory. There will be a file written in the .tsv file format for every single genome in the pangenome.
The columns of this file are described in the following table: 

| Column               | Description                                                                                                                  |
|----------------------|------------------------------------------------------------------------------------------------------------------------------|
| gene                 | the unique identifier of the gene                                                                                            |
| contig               | the contig that the gene is on                                                                                               |
| start                | the start position of the gene                                                                                               |
| stop                 | the stop position of the gene                                                                                                |
| strand               | The strand that the gene is on                                                                                               |
| ori                  | Will be T if the gene name is dnaA                                                                                           |
| family               | the family identifier to which the gene belongs to                                                                           |
| nb_copy_in_org       | The number of copy of the family in the organism (basically, if 1, the gene has no closely related paralog in that organism) |
| partition            | the partition to which the gene family of the gene belongs to                                                                |
| persistent_neighbors | The number of neighbors classified as 'persistent' in the pangenome graph                                                    |
| shell_neighbors      | The number of neighbors classified as 'shell' in the pangenome graph                                                         |
| cloud_neighbors      | The number of neighbors classidied as 'cloud' in the pangenome graph                                                         |

Those files can be generated as such: 

`ppanggolin write -p pangenome.h5 --projection`


#### Gene families to genes associtations

You can write a list containing the gene family assigned to every single gene of your pangenome, in a file format extactly like the one provided by [MMseqs2](https://github.com/soedinglab/MMseqs2) through its subcommand 'createtsv'.
It is basically a three-column file listing the gene family name in the first column, and the gene names in the second. A third column is either empty, or has an "F" in it. In that case it indicates that the gene is potentially a gene fragment and not complete. This will be indicated only if the [defragmentation](https://github.com/labgem/PPanGGOLiN/wiki/PPanGGOLiN---step-by-step-pangenome-analysis#defragmentation) pipeline is used.

You can obtain it as such:  

`ppanggolin write -p pangenome.h5 --families_tsv`

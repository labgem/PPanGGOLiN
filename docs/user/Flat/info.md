When computing a pangenome, all of the information about it is saved in the .h5 file, notably parameters used at each step and metrics about the pangenome. You can easily retrieve those informations using the 'info' module.

This command prints information on stdout, and does not write any file.

### Content

This option indicates the following metrics about your pangenome, if they have been computed:
* The total number of genes
* The number of genomes
* The number of gene families
* The number of edges in the pangenome graph
* The number of persistent genes, with the minimal, maximal, sd and mean presence thresholds of the families in this partition
* The number of shell genes, with the minimal, maximal, sd and mean presence thresholds of the families in this partition
* The number of cloud genes, with the minimal, maximal, sd and mean presence thresholds of the families in this partition
* The number of partitions

Additionally, if you have predicted RGPs and spots (with the subcommands 'panrgp', 'rgp' and 'spot', or 'all'), you will have the following metrics:
* The number of RGPs (Regions of Genomic Plasticity)
* The number of spots of insertion

Additionally, if you have predicted modules (with the subcommands 'panmodule', 'module' or 'all'):
* The number of modules
* The number of gene families in modules

It is used as such:
` ppanggolin info -p pangenome.h5 --content` 

### Parameters

This option indicates, for each steps of the analysis, the PPanGGOLiN parameters that were used and the source of the data if appropriate.

It is used as such:

`ppanggolin info -p pangenome.h5 --parameters`

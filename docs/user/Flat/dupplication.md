This file is a .tsv file, with a single parameter written as a comment at the beginning of the file, which indicates the proportion of genomes in which a gene family must be present more than once to be considered 'duplicated' (and not single copy marker).
This file lists the gene families, their duplication ratio, their mean presence in the pangenome and whether it is considered a 'single copy marker' or not, which is particularly useful when calculating the completeness recorded in the [organisms statistics](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#organisms-statistics) file described previously.

It can be generated using the 'write' subcommand as such : 

`ppanggolin write -p pangenome.h5 --stats`

This command will also generate the 'organisms_statistics.tsv' file.

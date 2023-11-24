## PPanGGOLiN complete workflow analyses

We tried to make PPanGGOLiN relatively easy to use by making this **'complete workflow'** subcommand. 
It runs a pangenome analysis whose exact steps will depend on the input files you provide it with.
In the end, you will have a partitioned pangenome graph with predicted **RGP, spots and modules**. 


[//]: # (### PPanGGOLiN: Pangenome analyses from list of annotated files)

The minimal subcommand only need your own annotations files (using .gff or .gbff/.gbk files) 
as long as they include the genomic dna sequences, such as the ones provided by prokka.
 
```
ppanggolin all --anno organism.gbff.list
```

It uses parameters that we found to be generally the best when working with species pangenomes.

The file **organism.gbff.list** is a tab-separated file with the following organisation :

1. The first column contains a unique organism name
2. The second column the path to the associated annotation file
3. Each line represents an organism

An example with 50 _Chlamydia trachomatis_ genomes can be found in the [testingDataset/ directory](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/).

[//]: # (### PPanGGOLiN: Pangenome analyses from list of fasta files)
You can also give PPanGGOLiN .fasta files, such as:

```
ppanggolin all --fasta organism.fasta.list
```

Again you must use a tab-separated file but this time with the following organisation:

1. The first column contains a unique organism name
2. The second column the path to the associated FASTA file
3. Circular contig identifiers are indicated in the following columns
4. Each line represents an organism

Same, an example can be found in the [testingDataset/ directory](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/).


Congratulation, you build (maybe your first) pangenome graph and partitioned it in 3 different partition: **persistent**, **shell** and **cloud** (look at ([Gautreau et al. 2020](https://doi.org/10.1371/journal.pcbi.1007732)) for more information about partition). 
You also detect **RGP, spots and modules** in your pangenome.

The results of the workflow is saved in the  **pangenome.h5** file, which is in the HDF-5 file format.
When you run an analysis using this file as input, the results of that analysis will be added to the file to supplement the data that are already stored in it. 
The idea behind this is that you can store and manipulate your pangenome with PPanGGOLiN by using this file only. It will keep all the information about what was done, all the parameters used, and all the pangenome's content.

[//]: # (TODO add link)
```{tip}
Many option are available to tune your analysis. Take a look here.
```

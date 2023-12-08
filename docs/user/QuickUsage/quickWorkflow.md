## PPanGGOLiN complete workflow analyses

We tried to make PPanGGOLiN relatively easy to use by making a **'complete workflow'** subcommand calls `all`. 
It runs a pangenome analysis whose exact steps will depend on the input files you provide it with.
In the end, you will have a partitioned pangenome graph with predicted **RGP, spots and modules**. 


[//]: # (### PPanGGOLiN: Pangenome analyses from list of annotated files)

The minimal subcommand only need your own annotations files (using `.gff` or `.gbff`/`.gbk` files) 
as long as they include the genomic dna sequences, such as the ones provided by Prokka or Bakta.
 
```bash
ppanggolin all --anno organism.gbff.list
```

It uses parameters that we found to be generally the best when working with species pangenomes.

The file **organism.gbff.list** is a tab-separated file with the following organisation :

1. The first column contains a unique organism name
2. The second column the path to the associated annotation file
3. Each line represents an organism

An example with 50 _Chlamydia trachomatis_ genomes can be found in the [testingDataset](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list) directory.

[//]: # (### PPanGGOLiN: Pangenome analyses from list of fasta files)
You can also give PPanGGOLiN `.fasta` files, such as:

```
ppanggolin all --fasta organism.fasta.list
```

Again you must use a tab-separated file but this time with the following organisation:

1. The first column contains a unique organism name
2. The second column the path to the associated FASTA file
3. Circular contig identifiers are indicated in the following columns
4. Each line represents an organism

Same, an example can be found in the [testingDataset](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.fasta.list) directory.

```{tip}
To easily download genomes of a species of interest from the NCBI assemblies refseq or genbank, you can use the CLI tools  [NCBI Genome Downloading Scripts](https://github.com/kblin/ncbi-genome-download) or the [genome updater](https://github.com/pirovc/genome_updater) bash script. 

For instance to download with genome updater the GTDB refseq genomes of Bradyrhizobium japonicum, you can run the folowing command 
```bash
genome_updater.sh -d "refseq"  -o "B_japonicum_genomes" -M "gtdb" -T "s__Bradyrhizobium japonicum"
```

```


```{tip}
Downloading genomes from NCBI assemblies refseq or genbank for a species of interest can be easily accomplished using CLI tools like [NCBI Genome Downloading Scripts](https://github.com/kblin/ncbi-genome-download) or the [genome updater](https://github.com/pirovc/genome_updater) script.

For example, to download GTDB refseq genomes of Bradyrhizobium japonicum using genome updater, you can execute the following command:
```bash
genome_updater.sh -d "refseq" -o "B_japonicum_genomes" -M "gtdb" -T "s__Bradyrhizobium japonicum"
```

```


After the completion of the `all` command, a pangenome graph has been successfully constructed and partitioned into three distinct paritions: **persistent**, **shell**, and **cloud** (for further details on partitions, refer to [Gautreau et al. 2020](https://doi.org/10.1371/journal.pcbi.1007732)). Additionally, **RGPs, spots, and modules** have been detected within your pangenome.

The results of the workflow is saved in the  **pangenome.h5** file, which is in the HDF-5 file format.
When you run an analysis using this file as input, the results of that analysis will be added to the file to supplement the data that are already stored in it. 
The idea behind this is that you can store and manipulate your pangenome with PPanGGOLiN by using this file only. It will keep all the information about what was done, all the parameters used, and all the pangenome's content.

[//]: # (TODO add link)
```{tip}
Many option are available to tune your analysis. Take a look here.
```

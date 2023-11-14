# How to prepare your data for PPanGGOLiN

To build and partition a pangenome, PPanGGOLiN need a set of either DNA sequences or provided genome annotations. In order to help you to start with PPanGGOLiN, you can follow this step to download some genomes from _Bradyrhizobium japonicum_. These genomes will be our base line all along the tutorial. If you already have your genome you can directly go to [the input file creation](#create-your-list-of-genomes-file)

## Get B. _japonicum_ genomics data

```{tip}
To download our genomes, we are going to use [genome_updater](https://github.com/pirovc/genome_updater).
Other solution exist such as [ncbi genome downloading scripts](https://github.com/kblin/ncbi-genome-download). Feel free to use the best and easiest way for you.
```

### GTDB genomes

To obtain the genomes of B. _japonicum_ from the [GTDB database](https://gtdb.ecogenomic.org/), you must use the name of the species in GTDB.

```
genome_updater.sh -d "refseq,genbank" -f "genomic.gbff.gz" -o "B_japonicum_genomes" -M "gtdb" -T "s__Bradyrhizobium japonicum"
```

### NCBI/GenBank genomes

To obtain the genomes of B. _japonicum_ from the [NCBI](https://www.ncbi.nlm.nih.gov/), you must use its taxonomic ID.

```
genome_updater.sh -d "refseq,genbank" -f "genomic.gbff.gz" -o "B_japonicum_genomes" -M "ncbi" -T "375"
```

## Create your list of genomes file

PPanGGOLiN use the list of genomes as input for some command, such as the workflow.
The file is a tsv-separated file with the following organisation :

1. The first column contains a unique organism name
2. The second column the path to the associated annotated file
3. Each line represents an organism

```{note}
It's also possible to use fasta file as input.
Look at the documentation.
```

If you are using the annotated genomes (*GBFF*, *GFF*, *GBK*), you can generate your file with the following command

```
for file in $(ls B_japonicum_genomes/*/files/*.gz);do genome=$(echo $file | cut -d'/' -f4 | cut -d'_' -f1-3); echo -e "$genome\t$file"; done > organism_gbff.list      
```

**You're now ready to build the pangenome !!!**
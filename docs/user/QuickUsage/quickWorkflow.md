## PPanGGOLiN complete workflow analyses

We tried to make PPanGGOLiN relatively easy to use by making a **'complete workflow'** subcommand called `all`. 
It runs a pangenome analysis whose exact steps will depend on the input files you provide it with.
In the end, you will have a partitioned pangenome graph with predicted **RGP, spots and modules**. 

[//]: # (### PPanGGOLiN: Pangenome analyses from list of annotated files)



```{mermaid}

---
title: "Workflow Overview: Steps launched by the all command"
align: center
---

%%{init: {'theme':'default'}}%%


graph LR

    i[input genomes] --> a
   

        r:::panrgp
        s:::panrgp
        m:::panmodule

        subgraph Pangenome creation
            a:::workflow
            c:::workflow
            g:::workflow
            p:::workflow
            a("annotate") --> c
            c(cluster) --> g(graph)
            g(graph) --> p(partition)
        end

        subgraph Functional module
        p --> m(module)
        end
        
        subgraph Region of Genomic Plasticity
        
        p --> r(rgp)
        r --> s(spot)
        end

    p --> f[pangenome.h5]
    s --> f
    m --> f

        
    classDef panrgp fill:#84d191
    classDef panmodule fill:#d44066
    classDef workflow fill:#d4ae40


```

The minimal subcommand only need your own annotations files (using `.gff` or `.gbff`/`.gbk` files) 
as long as they include the genomic dna sequences, such as the ones provided by Prokka or Bakta.
 
```bash
ppanggolin all --anno genomes.gbff.list
```

It uses parameters that we found to be generally the best when working with species pangenomes.

The file **genomes.gbff.list** is a tab-separated file with the following organisation :

1. The first column contains a unique genome name
2. The second column the path to the associated annotation file
3. Each line represents a genome

An example with 50 _Chlamydia trachomatis_ genomes can be found in the [testingDataset](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/genomes.gbff.list) directory.

[//]: # (### PPanGGOLiN: Pangenome analyses from list of fasta files)
You can also give PPanGGOLiN `.fasta` files, such as:

```
ppanggolin all --fasta genomes.fasta.list
```

Again you must use a tab-separated file but this time with the following organisation:

1. The first column contains a unique genome name
2. The second column the path to the associated FASTA file
3. Circular contig identifiers are indicated in the following columns
4. Each line represents a genome

Same, an example can be found in the [testingDataset](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/genomes.fasta.list) directory.

```{tip}
Downloading genomes from NCBI refseq or genbank for a species of interest can be easily accomplished using CLI tools like [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) or the [genome updater](https://github.com/pirovc/genome_updater) script.

For instance to download the GTDB refseq genomes of Bradyrhizobium japonicum with genome updater, you can run the following command 
```bash
genome_updater.sh -d "refseq"  -o "B_japonicum_genomes" -M "gtdb" -T "s__Bradyrhizobium japonicum"
```


After the completion of the `all` command, all of your genomes have had their genes predicted, the genes have been clustered into gene families, a pangenome graph has been successfully constructed and partitioned into three distinct paritions: **persistent**, **shell**, and **cloud**. Additionally, **RGP, spots, and modules** have been detected within your pangenome.

The results of the workflow is saved in the  **pangenome.h5** file, which is in the HDF-5 file format.
When you run an analysis using this file as input, the results of that analysis will be added to the file to supplement the data that are already stored in it. 
The idea behind this is that you can store and manipulate your pangenome with PPanGGOLiN by using this file only. It will keep all the information about what was done, all the parameters used, and all the pangenome's content.

```{tip}
Many option are available to tune your analysis. Take a look [here](../PangenomeAnalyses/pangenomeAnalyses.md#workflow).
```

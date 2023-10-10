(step-by-step-section)=
# Pangenome graph building and partitionning
The workflow subcommand of PPanGGOLiN uses what you can use through the other subcommands. 
In this section, the different steps will be described along with the inputs that can be given, and which subcommands to use to run those specific parts yourself.

This is useful only if you want to customize your workflow parameters. 
If you want to build the pangenome of your species without tuning parameters, you can use the subcommand `ppanggolin workflow` as described in the [introduction](#basic).

(annotation)=
## Genomes annotation and storage

The first PPanGGOLiN step consist to get all information about genomes and store them in the pangenome file. 
PPanGGOLiN will get information about: genes, RNA, contigs and genomes.

```{include} step-by-step/annotation.md
```

## Clustering genes in gene families

In order to build the pangenome graph, we need to compute clusters of similar genes called gene families. 

```{include} step-by-step/clustering.md
```

## Building the pangenome graph

Now that we get the gene families we can build the pangenome graph.

```{include} step-by-step/graph.md
```

## Pangenome graph partitioning

The final step is to partition the pre-computed pangenome graph.

```{include} step-by-step/partition.md
```
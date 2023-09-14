(step-by-step-section)=
# Pangenome graph building and partitionning
The workflow subcommand of PPanGGOLiN uses what you can use through the other subcommands. In this section, the different steps will be briefly described along with the inputs that can be given, and which subcommands to use to run those specific parts yourself.

This is useful only if you want to customize your workflow parameters. If you want to build the pangenome of your species without tuning parameters, you can use the subcommand `ppanggolin workflow` as described in the introduction.

(annotation)=
## Genomes annotation and storage
```{include} step-by-step/annotation.md
```

## Clustering genes in gene families
```{include} step-by-step/clustering.md
```

## Building the pangenome graph
```{include} step-by-step/graph.md
```

## Pangenome graph partitioning 
```{include} step-by-step/partition.md
```
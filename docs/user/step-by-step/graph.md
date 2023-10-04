
To partition a pangenome graph, you need to build a said pangenome graph. 
This can be done through the `graph` subcommand. 
This will take a pangenome .h5 file as input and compute edges to link gene families together based on the genomic neighborhood.
The graph is constructed using the following subcommand : 

```
ppanggolin graph -p pangenome.h5
```

This subcommand has only a single other option, which is `-r` or `--remove_high_copy_number`. 
If used, it will remove the gene families that are too duplicated in your genomes.
This is useful if you want to visualize your pangenome afterward and want to remove the biggest hubs to have a clearer view. 
It can also be used to limit the influence of very duplicated genes such as transposase or ABC transporters in the partition step.


The resulting pangenome graph is saved in the pangenome.h5 file given as input.
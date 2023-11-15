<!-- ### Genes table with pangenome annotations -->

The `--table` option generates a TSV file for each genome, providing pangenome annotations for the genes.  These files are stored within a directory named 'tables'.


The table below outlines the columns found in these generated files:

| Column               | Description                                                                |
|----------------------|----------------------------------------------------------------------------|
| gene                 | Unique identifier for the gene                                             |
| contig               | Contig on which the gene is located                                        |
| start                | Start position of the gene                                                 |
| stop                 | Stop position of the gene                                                  |
| strand               | Strand on which the gene is on                                       |
| family               | Id of the gene's associated family within the pangenome             |
| nb_copy_in_org       | Number of family copies present in the organism; 1 indicates no close paralogs |
| partition            | Partition to which the gene family belongs in the pangenome                  |
| persistent_neighbors | Number of neighbors classified as 'persistent' in the pangenome graph        |
| shell_neighbors      | Number of neighbors classified as 'shell' in the pangenome graph             |
| cloud_neighbors      | Number of neighbors classified as 'cloud' in the pangenome graph             |
| RGP                  | Name of the Region of Genomic Plasticity (RGP) if the gene is within an RGP  (present only if RGPs have been predicted) |
| spot                 | Spot ID in which the RGP is inserted (present only if RGPs and spot have been predicted)   |
| module               | Module ID of the gene family (present if modules have been predicted)             |


```{note}
Columns such as RGP, spot, and module are included only when these elements have been predicted in the pangenome.
```

Those files can be generated as such : 

```
ppanggolin write_genomes -p pangenome.h5 --table -o write_genomes_output
```


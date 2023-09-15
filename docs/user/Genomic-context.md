# Genomic context searching

From version 1.2.45, it is possible to search genomic context in a pangenome graph using PPanGGOLiN. A genomic context corresponds to a group of genes/proteins with a functional interest, often found together in the genomes. They are detected by extracting a subgraph obtained by filtering edges connecting the sequences of interest in the pangenome.

The analysis can be run on your formerly computed pangenomes and one or more genomic context. 

To search your genomic context of interest, there are two possibilities. You can search directly with genes/proteins sequences in a fasta file or use a list of the gene family ID. Both possibilities can be run in the same subcommand `context` and all the options are for tuning the parameters for the analysis.

## Search Genomic context with sequences

The first possibility to search your genomic context, you can use a fasta file with genes or proteins. The command can be launched as such:

`ppanggolin context -p pangenome.h5 --sequences protein.fasta`

This will search the genomic context in the computed pangenome and export the result in a tsv file.

To use this subcommand, be sure that your pangenome have gene families representatives sequences associated to it.

## Search with gene family ID.

The second possibility is to give a list of gene families ID used to compute the pangenome. You can run the subcommand like this:

`ppanggolin context -p pangenome.h5 --family families.txt`

This will search the common connected components in the computed pangenome and export the result in a tsv file.

In this case, you can give a pangenome without gene families representatives sequences. This option is compatible with a pangenome computed with an external clustering (see [cluster](https://github.com/labgem/PPanGGOLiN/wiki/PPanGGOLiN---step-by-step-pangenome-analysis#clustering) subcommand).

## Output format

In case of you are using families ID, you will only have as output the `gene_context.tsv` file. In the other case, you use sequences, you will have another output file to report the alignment between sequences and pangenome families (see detail in align subcommand).

There are 4 columns in `gene_context.tsv`. 

1. **geneContext ID**: identifier of the found context. It is incrementally generated, beginning with 1
2. **Gene family name**: Identifier of the gene family, from the pangenome, correspond to the found context
3. **Sequence ID**: Identifier of the searched sequence in the pangenome
4. **Nb Genomes**: Number of genomes where the genomic context is found
5. **Partition**: Partition of the gene family corresponding to the found context

In **sequence Id**, it is possible to find a NA value. This case, correspond to another gene family found in the context.

## Detailed options
| option name      | Description                                                                                                                                                                                                       |
|------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --no_defrag      | Do not use the defragmentation step, to align sequences with MMseqs2                                                                                                                                              |
| --identity       | Minimum identity percentage threshold                                                                                                                                                                             |
| --coverage       | Minimum coverage percentage threshold                                                                                                                                                                             |
| -t, --transitive | Size of the transitive closure used to build the graph. This indicates the number of non-related genes allowed in-between two related genes. Increasing it will improve precision but lower sensitivity a little. |
| -s, --jaccard    | Minimum jaccard similarity used to filter edges between gene families. Increasing it will improve precision but lower sensitivity a lot.                                                                          |
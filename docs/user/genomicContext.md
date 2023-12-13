# Prediction of Genomic Context

The PPanGGOLiN `context` command enables the identification of genomic contexts for query proteins. These contexts consist of genes commonly found in proximity to the proteins of interest in the different genomes.

The analysis can be run on a formerly computed pangenomes and users can query one or multiple genes at once. The search can be conducted either directly with gene/protein sequences in a FASTA file or by utilizing a list of gene family IDs. Both methods are seamlessly integrated within the `context` subcommand.


## Search Genomic context with sequences

To search your genomic context, you can use a FASTA file containing genes or proteins. The command can be launched as such:

`ppanggolin context -p pangenome.h5 --sequences protein.fasta`

To utilize this subcommand, ensure that your pangenome contains sequences associated with gene family representatives. This is the case with a pangenome computed with an external clustering (see the [cluster](./PangenomeAnalyses/pangenomeCluster.md) subcommand).

## Search with gene family ID.

Another possibility is to give a list of gene families ID used to compute the pangenome. You can run the subcommand like this:

`ppanggolin context -p pangenome.h5 --family families.txt`


In this scenario, you can give a pangenome without gene families representatives sequences. This option is compatible with a pangenome computed with an external clustering (see the [cluster](./PangenomeAnalyses/pangenomeCluster.md) subcommand).

## Output format

If you are using families IDs, the only output you will receive is the `gene_context.tsv` file. If you use sequences, you will have another output file that report the alignment between sequences and pangenome families (see detail in [align subcommand](align.md#align-external-genes-to-a-pangenome)).

There are 6 columns in `gene_context.tsv`. 

1. **geneContext_ID**: Identifier of the found context. It is incrementally generated, beginning with 1
2. **Gene_family_name**: Identifier of the gene family, from the pangenome, correspond to the found context
3. **Sequence_ID**: Identifier of the searched sequence in the pangenome
4. **Nb_Genomes**: Number of genomes where the genomic context is found
5. **Partition**: Partition of the gene family corresponding to the found context
6. **Target_family**: Whether the family is a target family, meaning it matches an input sequence, or a family provided as input.

```{note}
In **sequence_ID**, it is possible to find a NA value. This corresponds to cases where a gene family other than the one specified by the user is found in the context.
```


## Detailed options

| option name | Description |
|-----------------------------|---------------------------------------------------------------------------|
| --fast | Use representative sequences of gene families for input gene alignment. This option is recommended for faster processing but may be less sensitive. By default, all pangenome genes are used for alignment. This argument makes sense only when --sequence is provided. (default: False) |
| --no_defrag | Do not use the defragmentation step, to align sequences with MMseqs2 (default: False) |
| --identity | Minimum identity percentage threshold (default: 0.8)|
| --coverage | Minimum coverage percentage threshold (default: 0.8)|
| -t, --transitive | Size of the transitive closure used to build the graph. This indicates the number of non-related genes allowed in-between two related genes. Increasing it will improve precision but lower sensitivity a little. (default: 4) |
| -s, --jaccard | Minimum jaccard similarity used to filter edges between gene families. Increasing it will improve precision but lower sensitivity a lot. (default: 0.85) |
| -w, --window_size | Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. (default: 5) |

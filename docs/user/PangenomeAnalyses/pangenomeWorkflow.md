PPanGGOLiN was created to assist users in their analysis while ensuring ease of use. 
This was achieved by incorporating a 'workflow' command that enables the construction and partitioning of a pangenome using genomic data. 
This command launches the [annotation](./pangenomeAnalyses.md#annotation), [clustering](./pangenomeAnalyses.md#clustering), [graph](./pangenomeAnalyses.md#graph) and [partition](./pangenomeAnalyses.md#partition) commands described below.

You need to provide a tab-separated list of either annotated files or fasta files. The expected format is detailed [here](./pangenomeAnalyses.md#annotation)

You can use the workflow with annotated files such as: 
```
ppanggolin workflow --anno organism.gbff.list
```

For fasta files, you have to change for this command: 
```
ppanggolin workflow --fasta organism.fasta.list
```

Moreover, as explained [here](./pangenomeAnalyses.md#read-clustering), it is also possible to provide your own clustering in the workflow command as such:

```
ppanggolin workflow --anno organism.gbff.list --clusters clusters.tsv
```

All options are common to the command explain below, so you can have complete information by reading them.
The only workflow specific option is `--no_flat_files`. 
This option prevents the automatic generation of the output files listed and described [here](./pangenomeAnalyses.md#pangenome-outputs).
If you are not familiar with the outputs available in PPanGGOLiN, we recommend reading this section and 
 not using this option to have all the result automatically generated despite it cool be time consuming.

```{tip}
In the workflow CLI, it is not possible to tune all the options available in all the steps. 
For a fully optimized analysis, we recommend to use the configuration file as described [here](../practicalInformation.md#configuration-file)
```

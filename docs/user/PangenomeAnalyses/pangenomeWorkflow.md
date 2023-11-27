PPanGGOLiN was designed to help the users in their analyses and to be easy to use. 
For that we made a 'workflow' command to construct and partitioned a pangenome from genomic data. 
This command will launch the [annotation](./pangenomeAnalyses.md#annotation), [clustering](./pangenomeAnalyses.md#clustering), [graph](./pangenomeAnalyses.md#graph) and [partition](./pangenomeAnalyses.md#partition) command detailed below.

You need to provide a tab-separated list of either annotated files or fasta files. The expected format is detailed [here](./pangenomeAnalyses.md#annotation)

You can use the workflow with annotated files such as: 
```
ppanggolin workflow --anno organism.gbff.list
```

For fasta files, you have to change for this command: 
```
ppanggolin workflow --fasta organism.fasta.list
```

Moreover, as its explained [here](./pangenomeAnalyses.md#read-clustering), it's also possible to provide your own clustering in the workflow command as such:

```
ppanggolin workflow --anno organism.gbff.list --clusters clusters.tsv
```

All the option are common to the command explain below, so you can have complete information by reading them. 
There is only one option workflow specific: `--no_flat_files`. 
This option will avoid to automatically generate the output files listed and described [here](./pangenomeAnalyses.md#pangenome-outputs).
If you're not familiar with the outputs available in PPanGGOLiN, we recommend to read this part and 
to not use this option to have all the result automatically generated despite it cool be time consuming.

```{tip}
In the workflow CLI, it's not possible to tune all the option available in all the step. 
If you want to have a complete tuned analyses we recommend to use the configuration file as describe [here](../practicalInformation.md#configuration-file)
```
PPanGGOLiN was created with the idea to make it both easy to use for beginners, and fully customizable for experts.
Ease of use has been achieved by incorporating a workflow command that allows the construction and partitioning of a pangenome using genomic data.
The command has only one mandatory option, and predefined parameters adapted to pangenomes at the scale of a bacterial species.
This command launches the [annotation](./pangenomeAnalyses.md#annotation), [clustering](./pangenomeCluster.md#cluster-genes-into-gene-families), [graph](./pangenomeAnalyses.md#graph) and [partition](./pangenomeAnalyses.md#partition) commands described below.

<br>
<br>

```{mermaid}

---
title: "Workflow Overview: Steps launched by the workflow command"
align: center
---

%%{init: {'theme':'default'}}%%


graph LR

    i[input genomes] --> a


        subgraph Pangenome creation
            a:::workflow
            c:::workflow
            g:::workflow
            p:::workflow
            a("annotate") --> c
            c(cluster) --> g(graph)
            g(graph) --> p(partition)
        end

    p --> f[pangenome.f5]

        
    classDef panrgp fill:#4066d4
    classDef panmodule fill:#d44066
    classDef workflow fill:#d4ae40


```

<br>
<br>

To use this command, you need to provide a tab-separated list of either annotation files (gff3 or gbff) or fasta files. The expected format is detailed [in the annotation section](./pangenomeAnalyses.md#annotation)

You can use the workflow with annotation files as such: 
```
ppanggolin workflow --anno organism.gbff.list
```

For fasta files, you have to change for: 
```
ppanggolin workflow --fasta organism.fasta.list
```

Moreover, as detailed [in the section about providing your gene families](./pangenomeAnalyses.md#read-clustering), 
if you wish to use different gene clustering methods than those provided by PPanGGOLiN,
it is also possible to provide your own clustering results with the workflow command as such:

```
ppanggolin workflow --anno organism.gbff.list --clusters clusters.tsv
```

All the workflow parameters are obtained from the commands explained below, except for the `--no_flat_files` option, which solely pertains to it. This option prevents the automatic generation of the output files listed and described [in the pangenome output section](./pangenomeAnalyses.md#pangenome-outputs).
If you are unfamiliar with the output available in PPanGGOLiN, we recommend that you do not use this option, so that all results are automatically generated, even though this may take some time.

```{tip}
In the workflow CLI, it is not possible to tune all the options available in all the steps. 
For a fully optimized analysis, you can either launch the subcommands one by one as described below, or you can use the configuration file as described [here](../practicalInformation.md#configuration-file)
```





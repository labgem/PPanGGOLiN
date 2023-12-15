## Purpose

Regions of Genome Plasticity (RGPs) are clusters of genes made of shell and cloud genomes in the pangenome graph. Most of them arise from Horizontal gene transfer (HGT) and correspond to Genomic Islands (GIs). RGPs from different genomes can be grouped in spots of insertion based on their conserved flanking persistent genes, rather than their gene content, to find out which are located in the same locations in the genome. The panRGP methods and its subcommands and subsequent output files are made to detect describe as thoroughly as possible those Regions of Genome Plasticity across all genomes of the pangenome.

Those methods were supported by the [panRGP publication](https://doi.org/10.1093/bioinformatics/btaa792) which can be read to have their methodological descriptions and justifications.


## PanRGP

This command works exactly like [workflow](../PangenomeAnalyses/pangenomeAnalyses.md#workflow). The difference is that it will run additional analyses to characterize Regions of Genome Plasticity (RGP).



```{mermaid}

---
title: "Workflow Overview: Steps launched by the panrgp command"
align: center
---

%%{init: {'theme':'default'}}%%


graph LR

    i[input genomes] --> a
   

        r:::panrgp
        s:::panrgp

        subgraph Pangenome creation
            a:::workflow
            c:::workflow
            g:::workflow
            p:::workflow
            a("annotate") --> c
            c(cluster) --> g(graph)
            g(graph) --> p(partition)
        end
        
        subgraph Region of Genomic Plasticity
        
        p --> r(rgp)
        r --> s(spot)
        end

    p --> f[pangenome.f5]
    s --> f[pangenome.f5]

        
    classDef panrgp fill:#84d191
    classDef panmodule fill:#d44066
    classDef workflow fill:#d4ae40


```


You can use the `panrgp` with annotation (gff3 or gbff) files with `--anno` option, as such: 
```bash
ppanggolin panrgp --anno organism.gbff.list
```

For fasta files, you need to use the alternative `--fasta` option, as such:
```bash
ppanggolin panrgp --fasta organism.fasta.list
```

Just like [workflow](../PangenomeAnalyses/pangenomeAnalyses.md#workflow), this command will deal with the [annotation](../PangenomeAnalyses/pangenomeAnalyses.md#annotation), [clustering](../PangenomeAnalyses/pangenomeAnalyses.md#compute-pangenome-gene-families), [graph](../PangenomeAnalyses/pangenomeAnalyses.md#graph) and [partition](../PangenomeAnalyses/pangenomeAnalyses.md#partition) commands by itself.
Then, the RGP detection is ran using [rgp](#rgp-detection) after the pangenome partitionning. Once all RGP have been computed, those found in similar genomic contexts in the genomes are gathered into spots of insertion using [spot](#spot-prediction).

If you want to tune the rgp detection, you can use the `rgp` command after the `workflow` command. If you wish to tune the spot detection, you can use the `spot` command after the `rgp` command. Additionally, you have the option to utilize a configuration file to customize each detection within the `panrgp` command.

More detail about RGP detection and the spot of insertion prediction can be found in the [panRGP publication](https://doi.org/10.1093/bioinformatics/btaa792)

## RGP detection

Once partitions have been computed, you can predict the Regions of Genome Plasticity. 
This subcommand's options are about tuning parameters for the analysis. 

You can do it as such:

```bash
ppanggolin rgp -p pangenome.h5
```

This will predict RGP and store results in the HDF5 file.

This command has 5 options that will directly impact the method's results. If you are not sure to understand what they do, you should not change their values as they'll probably behave exactly the way you want them to.
The options `--variable_gain`, `--persistent_penalty` are both linked to each other. They will impact respectively by how much a persistent gene will penalize the prediction of a RGP, and oppositely how a variable gene (shell, cloud, or multigenic persistent) will promote the prediction of a RGP. The `--min_score` is the filter that will classify a genome region as RGP or not.
Users looking to change those 3 parameters should consider reading the Materials and Methods of the [panRGP publication](https://doi.org/10.1093/bioinformatics/btaa792) to fully understand how those values impact the prediction.

The two other options are more straightforward. The `--min_length` will indicate the minimal size in base pair that a RGP should be to be predicted. The `--dup_margin` is a filter used to identify persistent gene families to consider as multigenic. Gene families that have more than one gene in more than `--dup_margin` genomes will be classified as multigenic, and as such considered as "variable" genes.

After this command is executed, a single output file that will list all of the predictions can be written, the [regions_of_genomic_plasticity.tsv](./rgpOutputs.md#rgp) file.

## Spot prediction

To study RGP that are found in the same area in different genomes, we gather them into 'spots of insertion'. These spots are groups of RGP that do not necessarily have the same gene content but have similar bordering _persistent_ genes. We run those analyses to study the dynamic of gene turnover of large regions in bacterial genomes. In this way, spots of the same pangenome can be compared and their dynamic can be established by comparing their different metrics together. Detailed descriptions of these metrics can be found in the [RGP and spot output section](./rgpOutputs.md#rgp-outputs).

Spots can be computed once RGP have been predicted. You can do that using:

```bash
ppanggolin spot -p pangenome.h5
```

This command has 3 options that can change its results:

- `--set_size` defines the number of the bordering _persistent_ genes that will be compared. (default: 3)
- `--overlapping_match` defines the minimal number of bordering persistent genes that must be identical and in the same order to consider two RGP borders as being similar. This means that, when comparing two borders, if they overlap by this amount of persistent genes they will be considered similar, even if they start or stop with different genes.  (default: 2)
- `--exact_match_size` Defines the minimal number of firstly bordering persistent genes to consider two RGP borders as being similar. If two RGP borders start with this amount of identical persistent gene families, they are considered similar. (default: 1)

The two other options are related to the 'spot graph' used to predict spots of insertion.

- `--spot_graph` writes the spot graph once predictions are computed
- `--graph_formats` defines the format in which the spot graph is written. Can be either gexf or graphml. (default: gexf)

You can the use the dedicated subcommand [draw](./rgpOutputs.md#draw-spots) to draw interactive figures for any given spot with the python library [bokeh](http://docs.bokeh.org/en/latest/). Those figures can can be visualized and modified directly in the browser. This plot is described [here](./rgpOutputs.md#draw-spots)

Multiple files can then be written describing the predicted spots and their linked RGP, such as a [file linking RGPs with their spots](./rgpOutputs.md#spots) and a [file showing multiple metrics for each spot](./rgpOutputs.md#summarize-spots).

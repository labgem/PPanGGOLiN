PPanGGOLiN provides multiple outputs to describe a pangenome. In this section the different outputs will be described.

In most cases it will provide with a HDF-5 file named "pangenome.h5". This file stores all the information about your pangenome and the analysis that were run. If given to ppanggolin through most of the subcommands, it will read information from it. This is practical as you can regenerate figures or output files, or rerun parts of the analysis without redoing everything.

In this section, each parts will describe a possible output of PPanGGOLiN, and will be commented with the command line that generates it using the HDF5 file, which is assumed to be called 'pangenome.h5'.

When using the same subcommand (like 'write' or 'draw' that can help you generate multiple file each), you can provide multiple options to write all of the file formats that you desire at once.

# Draw

## U-shaped plot

A U-shaped plot is a figure presenting the number of families (y axis) per number of organisms (x axis).
It is a .html file that can be opened with any browser and with which you can interact, zoom, move around, mouseover to see numbers in more detail, and you can save what you are seeing as a .png image file.

It can be generated using the 'draw' subcommand as such : 

`ppanggolin draw -p pangenome.h5 --ucurve`

## tile plot

A tile plot is a heatmap representing the gene families (y axis) in the organisms (x axis) making up your pangenome. The tiles on the graph will be colored if the gene family is present in an organism and uncolored if absent. The gene families are ordered by partition, and the genomes are ordered by a hierarchical clustering based on their shared gene families (basically two genomes that are close together in terms of gene family composition will be close together on the figure).

This plot is quite helpful to observe potential structures in your pangenome, and can also help you to identify eventual outliers. You can interact with it, and mousing over a tile in the plot will indicate to you which is the gene identifier(s), the gene family and the organism that corresponds to the tile.

If you build your pangenome using the 'workflow' subcommand and you have more than 500 organisms, only the 'shell' and the 'persistent' partitions will be drawn, leaving out the 'cloud' as the figure tends to be too heavy for a browser to open it otherwise.

It can be generated using the 'draw' subcommand as such : 

`ppanggolin draw -p pangenome.h5 --tile_plot`

and if you do not want the 'cloud' gene families as it is a lot of data and can be hard to open with a browser sometimes, you can use the following option : 

`ppanggolin draw -p pangenome.h5 --tile_plot --nocloud`


## Spot plots

For versions 1.2.30 and above, the 'draw' command can draw specific spots of interest, whose ID are provided, or all the spots if you wish.
It will also write a gexf file, which corresponds to the gene families and their organization within the spots. It is basically a subgraph of the pangenome, consisting of the spot itself.
The command can be used as such:

`ppanggolin draw -p pangenome.h5 --spots all` will draw an interactive .html figure and a gexf file for all the spots.

If you are interested in only a single spot, you can use its identifier to draw it, as such:

`ppanggolin draw -p pangenome.h5 --spots spot_34` for spot_34, for example.

The interactive figures that are drawn look like this:

![interactive figure](https://github.com/labgem/PPanGGOLiN/raw/master/images/drawspot_example.png)

The plot represents the different gene organizations that are found in the spot. If there are RGPs with identical gene organization, the organization is represented only once (the represented RGP is picked at random among all identical RGPs). The list of RGPs with the same organization is accessible in the file written alongside the figure called 'spot_X_identical_rgps.tsv', with X the spot_id.

They can be edited using the sliders and the radio buttons, to change various graphical parameters, and then the plot itself can be saved using the save button on the right of the screen, if need be.
For the gexf file, you can see how to visualize it in the section about the [pangenome gexf](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#gexf-and-light-gexf).

# Rarefaction

This figure is not drawn by default in the 'workflow' subcommand as it requires a lot of computations. It represents the evolution of the number of gene families for each partition as you add more genomes to the pangenome. It has been used a lot in the literature as an indicator of the diversity that you are missing with your dataset on your taxonomic group. The idea is that if at some point when you keep adding genomes to your pangenome you do not add any more gene families, you might have access to your entire taxonomic group's diversity. On the contrary if you are still adding a lot of genes you may be still missing a lot of gene families. 

There are 8 partitions represented. For each of the partitions there are multiple representations of the observed data. You can find the observed means, medians, 1st and 3rd quartiles of the number of gene families per number of genome used. And you can find the fitting of the data by the Heaps' law, which is usually used to represent this evolution of the diversity in terms of gene families in each of the partitions.

It can be generated using the 'rarefaction' subcommand, which is dedicated to drawing this graph, as such : 

`ppanggolin rarefaction -p pangenome.h5`

A lot of options can be used with this subcommand to tune your rarefaction curves, most of them are the same as with the `partition` workflow.
The following 3 are related to the rarefaction alone:

- `--depth` defines the number of sampling for each number of organism (default 30)
- `--min` defines the minimal number of organisms in a sample (default 1)
- `--max` defines the maximal number of organisms in a sample (default 100)

So for example the following command:
`ppanggolin rarefaction -p pangenome.h5 --min 5 --max 50 --depth 30`

Will draw a rarefaction curve with sample sizes between 5 and 50 (between 5 and 50 genomes will be used), and with 30 samples at each point (so 30 samples of 5 genomes, 30 samples or 6 genomes ... up to 50 genomes).

# Write
## Organisms statistics

The organisms_statistics.tsv file is a tab-separated file describing the content of each of the genome used for building the pangenome. It might be useful when working with fragmented data such as MAGs or if you suspect some of your genomes to be chimeric, or to not belong to your taxonomic group (as those genomes will be outliers regarding to the numbers in this file).
The first lines starting with a '#' are indicators of parameters used when generating the numbers describing each organisms, and should not be read if loading this into a spreadsheet. They will be skipped automatically if you load this file with R.

This file is made of 15 columns described in the following table

| Column                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| organism               | Indicates the organism's name to whom the provided genome belongs to                                                                                                                                                                                                                                                                                                                                                                                                                    |
| nb_families            | Indicates the number of gene families present in that genome                                                                                                                                                                                                                                                                                                                                                                                                                            |
| nb_persistent_families | The number of persistent families present in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                |
| nb_shell_families      | The number of shell families present in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| nb_cloud_families      | The number of cloud families present in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| nb_exact_core          | The number of exact core families present in that genome. This number should be identical in all genomes.                                                                                                                                                                                                                                                                                                                                                                               |
| nb_soft_core           | The number of soft core families present in that genome. The threshold used is indicated in the #soft_core line at the beginning of the file, and is 0.95 by default.                                                                                                                                                                                                                                                                                                                   |
| nb_genes               | The number of genes in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| nb_persistent_genes    | The number of genes whose family is persistent in that genome                                                                                                                                                                                                                                                                                                                                                                                                                           |
| nb_shell_genes         | The number of genes whose family is shell in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                |
| nb_cloud_genes         | The number of genes whose family is cloud in that genome                                                                                                                                                                                                                                                                                                                                                                                                                                |
| nb_exact_core_genes    | The number of genes whose family is exact core in that genome                                                                                                                                                                                                                                                                                                                                                                                                                           |
| nb_soft_core_genes     | The number of genes whose family is soft core in that genome                                                                                                                                                                                                                                                                                                                                                                                                                            |
| completeness           | This is an indicator of the proportion of single copy markers in the persistent that are present in the genome. While it is expected to be relatively close to 100 when working with isolates, it may be particularly interesting when working with very fragmented genomes as this provide a *de novo* estimation of the completeness based on the expectation that single copy markers within the persistent should be mostly present in all individuals of the studied taxonomic group |
| nb_single_copy_markers | This indicates the number of present single copy markers in the genomes. They are computed using the parameter duplication_margin indicated at the beginning of the file. They correspond to all of the persistent gene families that are not present in more than one copy in 5% (or more) of the genomes by default.                                                                                                                                                                  |

It can be generated using the 'write' subcommand as such : 

`ppanggolin write -p pangenome.h5 --stats`

This command will also generate the 'mean_persistent_duplication.tsv' file.

## pangenomeGraph files
The pangenome's graph can be given through multiple data formats, in order to manipulate it with other softwares.

### gexf and light gexf
The Graph can be given through the .gexf and through the _light.gexf files. The _light.gexf file will contain the gene families as nodes and the edges between gene families describing their relationship, and the .gexf file will contain the same thing, but also include more informations about each gene and each relation between gene families. 
We have made two different files representing the same graph because, while the non-light file is exhaustive, it can be very heavy to manipulate and most of the information in it are not of interest to everyone. The _light.gexf file should be the one you use to manipulate the pangenome graph most of the time.

They can be manipulated and visualised through a software called [Gephi](https://gephi.org/), with which we have made extensive testings, or potentially any other softwares or libraries that can read gexf files such as [networkx](https://networkx.github.io/documentation/stable/index.html) or [gexf-js](https://github.com/raphv/gexf-js) among others. 

Using Gephi, the layout can be tuned as illustrated below:

![Gephi layout](https://github.com/labgem/PPanGGOLiN/raw/master/images/gephi.gif)

We advise the Gephi "Force Atlas 2" algorithm to compute the graph layout with "Stronger Gravity: on" and "scaling: 4000" but don't hesitate to tinker the layout parameters.

In the _light.gexf file : 
The nodes will contain the number of genes belonging to the gene family, the most commun gene name (if you provided annotations), the most common product name(if you provided annotations), the partitions it belongs to, its average and median size in nucleotids, and the number of organisms that have this gene family.

The edges contain the number of times they are present in the pangenome.

The .gexf non-light file will contain in addition to this all the information about genes belonging to each gene families, their names, their product string, their sizes and all the information about the neighborhood relationships of each pair of genes described through the edges.

The light gexf can be generated using the 'write' subcommand as such : 

`ppanggolin write -p pangenome.h5 --light_gexf`

while the gexf file can be generated as such : 

`ppanggolin write -p pangenome.h5 --gexf`

### json
The json's file content corresponds to the .gexf file content, but in json rather than gexf file format. It follows the 'node-link' format as shown in [this example](https://observablehq.com/@d3/force-directed-graph) in javascript, or as used in the [networkx](https://networkx.github.io/documentation/stable/reference/readwrite/json_graph.html) python library and it should be usable with both [D3js](https://d3js.org/) and [networkx](https://networkx.github.io/documentation/stable/index.html), or any other software or library that supports this format.

The json can be generated using the 'write' subcommand as such : 

`ppanggolin write -p pangenome.h5 --json`

## gene presence absence

This file is basically a presence absence matrix. The columns are the genomes used to build the pangenome, the lines are the gene families. The identifier of the gene family is the gene identifier chosen as a representative.
 There is a 1 if the gene family is present in a genome, and 0 otherwise. It follows the exact same format than the 'gene_presence_absence.Rtab' file that you get from the pangenome analysis software [roary](https://sanger-pathogens.github.io/Roary/)

It can be generated using the 'write' subcommand as such : 

`ppanggolin write -p pangenome.h5 --Rtab`

## matrix

This file is a .csv file following a format alike the  gene_presence_absence.csv file generated by [roary](https://sanger-pathogens.github.io/Roary/), and works with [scoary](https://github.com/AdmiralenOla/Scoary) if you want to do pangenome-wide association studies.

It can be generated using the 'write' subcommand as such : 

`ppanggolin write -p pangenome.h5 --csv`

## mean persistent duplication

This file is a .tsv file, with a single parameter written as a comment at the beginning of the file, which indicates the proportion of genomes in which a gene family must be present more than once to be considered 'duplicated' (and not single copy marker).
This file lists the gene families, their duplication ratio, their mean presence in the pangenome and whether it is considered a 'single copy marker' or not, which is particularly useful when calculating the completeness recorded in the [organisms statistics](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#organisms-statistics) file described previously.

It can be generated using the 'write' subcommand as such : 

`ppanggolin write -p pangenome.h5 --stats`

This command will also generate the 'organisms_statistics.tsv' file.

## partitions

Those files will be stored in the 'partitions' directory and will be named after the partition that they represent (like persistent.txt for the persistent partition). In each of those file there will be a list of gene family identifiers that correspond to the gene families belonging to that partition, one family per line, should you need it for your pipelines or during your analysis.

You can generate those files as such :  

` ppanggolin write -p pangenome.h5 --partitions`

## projection

This option writes in a 'projection' directory. There will be a file written in the .tsv file format for every single genome in the pangenome.
The columns of this file are described in the following table : 

| Column               | Description                                                                                                                    |
|----------------------|--------------------------------------------------------------------------------------------------------------------------------|
| gene                 | the unique identifier of the gene                                                                                              |
| contig               | the contig that the gene is on                                                                                                 |
| start                | the start position of the gene                                                                                                 |
| stop                 | the stop position of the gene                                                                                                  |
| strand               | The strand that the gene is on                                                                                                 |
| ori                  | Will be T if the gene name is dnaA                                                                                                              |
| family               | the family identifier to which the gene belongs to                                                                             |
| nb_copy_in_org       | The number of copy of the family in the organism (basically, if 1, the gene has no closely related paralog in that organism) |
| partition            | the partition to which the gene family of the gene belongs to                                                                  |
| persistent_neighbors | The number of neighbors classified as 'persistent' in the pangenome graph                                                      |
| shell_neighbors      | The number of neighbors classified as 'shell' in the pangenome graph                                                           |
| cloud_neighbors      | The number of neighbors classidied as 'cloud' in the pangenome graph                                                           |

Those files can be generated as such : 

`ppanggolin write -p pangenome.h5 --projection`

## Gene families and genes

You can write a list containing the gene family assigned to every single gene of your pangenome, in a file format extactly like the one provided by [MMseqs2](https://github.com/soedinglab/MMseqs2) through its subcommand 'createtsv'.
It is basically a three-column file listing the gene family name in the first column, and the gene names in the second. A third column is either empty, or has an "F" in it. In that case it indicates that the gene is potentially a gene fragment and not complete. This will be indicated only if the [defragmentation](https://github.com/labgem/PPanGGOLiN/wiki/PPanGGOLiN---step-by-step-pangenome-analysis#defragmentation) pipeline is used.

You can obtain it as such :  

`ppanggolin write -p pangenome.h5 --families_tsv`

## Plastic regions

This file is a tsv file that lists all of the detected Regions of Genome Plasticity. This requires to have run the RGP detection analysis by either using the `panrgp` command or the `rgp` command.

It can be written with the following command:
`ppanggolin write -p pangenome.h5 --regions`

The file has the following format :

| column | description |
|--------|-------------|
| region | a unique identifier for the region. This is usually built from the contig it is on, with a number after it|
|organism| the organism it is in. This is the organism name provided by the user.|
|start| the start position of the RGP in the contig|
|stop| the stop position of the RGP in the contig|
|genes| the number of genes included in the RGP|
|contigBorder| this is a boolean column. If the RGP is on a contig border it will be True, otherwise, it will be False. This often can indicate that, if an RGP is on a contig border it is probably not complete.|
|wholeContig| this is a boolean column. If the RGP is an entire contig, it will be True, and False otherwise. If an RGP is an entire contig it can possibly be a plasmid, a region flanked with repeat sequences or a contaminant|

## Spots

This is a tsv file with two column. It links the spots of 'summarize_spots' with the RGPs of 'plastic_regions'.

It is written with the following command:
`ppanggolin write -p pangenome.h5 --spots`

|column|description|
|------|------------|
|spot_id| The spot identifier (found in the 'spot' column of 'summarize_spots')|
|rgp_id| the RGP identifier (found in 'region' column of 'plastic_regions')|

## Summarize spots

This is a tsv file that will associate each spot with multiple metrics that can indicate the dynamic of the spot.

It is written with the following command:
`ppanggolin write -p pangenome.h5 --spots`

|column| description|
|-------|------------|
|spot| the spot identifier. It is unique in the pangenome|
|nb_rgp| the number of RGPs present in the spot|
|nb_families| The number of different gene families that are found in the spot|
|nb_unique_family_sets| The number of RGPs with different gene family content. If two RGPs are identical, they will be counted only once. The difference between this number and the one provided in 'nb_rgp' can be a strong indicator on whether their is a high turnover in gene content in this area or not|
|mean_nb_genes| the mean number of genes on RGPs in the spot|
|stdev_nb_genes| the standard deviation of the number of genes in the spot|
|max_nb_genes| the longest RGP in number of genes of the spot|
|min_nb_genes| the shortest RGP in number of genes of the spot|

## Borders

Each spot has at least one set of gene families bordering them. To write the list of gene families bordering a spot, you need to use the following option:
`ppanggolin write -p pangenome.h5 --borders`

It will write a .tsv file with 4 columns:

|column| description|
|-------|------------|
|spot_id| the spot identifier. It is unique in the pangenome|
|number| the number of RGPs present in the spot that have those bordering genes|
|border1| Comma-separated list of gene families of the 1st border|
|border2| Comma-separated list of gene families of the 2nd border|

As there can be some variation in the borders, some spots will have multiple borders and as such multiple lines in this file.
The sum of the number for each spot_id should be exactly the number of RGPs in the spot.

## Modules
### Functional modules
This .tsv file lists the modules and the gene families that belong to them. It lists one family per line, and there are multiple line for each module.
It is written along with other files with the following command:
`ppanggolin write -p pangenome.h5 --modules`

It follows the following format:
|column|description|
|------|------------|
|module_id| The module identifier|
|family_id| the family identifier|

### Modules in organisms
This .tsv file lists for each organism the modules that are present and how complete they are. Since there are some variability that are allowed in the module predictions, occasionnally some modules can be incomplete in some of the organisms where they are found.
This file is written along with other files with the following command:
`ppanggolin write -p pangenome.h5 --modules`

And it follows the following format:
|column|description|
|------|------------|
|module_id| The module identifier|
|organism| the organism which has the indicated module|
|completion| a value between 0.0 and 1.0 which indicates how complete (in terms of gene family) the module is in the given organism|

### modules summary
This .tsv file lists a few characteristics for each detected module. There is one line for each module.
The file is written along with other files with the following command:
`ppanggolin write -p pangenome.h5 --modules`

And it follows the following format:
|column|description|
|------|------------|
|module_id| The module identifier|
|nb_families| The number of families which are included in the module The families themselves are listed in the 'functional_modules.tsv' file.|
|nb_organisms|The number of organisms in which the module is found. Those organisms are listed in the 'modules_in_organisms.tsv' file.|
|partition| The average partition of the families in the module.|
|mean_number_of_occurrence| the mean number of time a module is present in each organism. The expected value is around one, but it can be more if it is a module often repeated in the genomes (like a phage).|

## spot modules
This command is available only if both modules and spots have been computed for your pangenome (see the command `all`, or the commands `spot` and `module` for that).
It indicates which modules are present in which spot and in which RGP.
The files are written with the following command:
```ppanggolin write -p pangenome.h5 --spot_modules```
The format of the 'modules_spots.tsv' file is the following:

|column|description|
|------|------------|
|module_id| The module identifier|
|spot_id| the spot identifier|

The file 'modules_RGP_lists.tsv' lists RGPs that have the same modules. Those RGPs can have different gene families, however they will not have any other module than those that are indicated. The format of the 'modules_RGP_lists.tsv' is the following:

|column|description|
|------|------------|
|representative_RGP| an RGP deemed representative for the group, and serving as a 'group of rgp id'(randomly picked)|
|nb_spots| The number of spots in which we see the RGPs which have the modules listed afterwards|
|mod_list| a list of the modules that are in the indicated RGPs|
|RGP_list| a list of RGP that include exactly the modules listed previously|

This information can also be visualized through figures that can be drawn with `ppanggolin draw --spots` (see [Spot plots](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#spot-plots), and which can display modules.

# Fasta
This command is available from 1.1.98 and on.
This command can be used to write fasta sequences of the pangenome or of parts of the pangenome. Most options expect a partition to write. Available partitions are:
* 'all' for the entire pangenome.
* 'Persistent' for persistent families
* 'Shell' for shell genes or families
* 'Cloud' for cloud genes or families
* 'rgp' for genes or families found in RGPs
* 'core' for core genes or families
* 'softcore' for softcore genes or families

When using the 'softcore' filter, the '--soft_core' option can be used to modily the threshold used to determine what is part of the soft core. It is at 0.95 by default.

## Genes

This option can be used to write the nucleotide CDS sequences. It can be used as such, to write all of the genes of the pangenome for example:

```ppanggolin fasta -p pangenome.h5 --output MY_GENES --genes all```

Or to write only the persistent genes:

```ppanggolin fasta -p pangenome.h5 --output MY_GENES --genes persistent```


## Protein families

This option can be used to write the protein sequences of the representative sequences for each family. It can be used as such for all families:

```ppanggolin fasta -p pangenome.h5 --output MY_PROT --prot_families all```

or for all of the shell families for example:

```ppanggolin fasta -p pangenome.h5 --output MY_PROT --prot_families shell```


## Gene families

This option can be used to write the gene sequences of the representative sequences for each family. It can be used as such:

```ppanggolin fasta -p pangenome.h5 --output MY_GENES --gene_families all```

or for the cloud families for example:

```ppanggolin fasta -p pangenome.h5 --output MY_GENES --gene_families cloud```

## Regions

This option can be used to write the nucleotide sequences of the detected RGPs. It requires the fasta sequences that were used to compute the pangenome as they were provided originally when you computed your pangenome. This command only has two filters:
* all, for all regions
* complete, for only the 'complete' regions which are not on a contig border

It can be used as such:

```ppanggolin fasta -p pangenome.h5 --output MYREGIONS --regions all --fasta organisms.fasta.list```

# MSA

This command is available from 1.1.103 and on.
It is used to call [mafft](https://mafft.cbrc.jp/alignment/software/) with default options to compute MSA of any partition of the pangenome. Using multiple cpus is recommended as it is quite demanding in computational resources.

By default it will write the strict 'core' (genes that are present in absolutely all genomes) and remove any duplicated genes. Beware however that, if you have many genomes (over 1000), the core will likely be either very small or even empty if you have fragmented genomes.

It will write one MSA for each family. You can then provide the directory where the MSA are written to [IQ-TREE](https://github.com/Cibiv/IQ-TREE) for example, to do phylogenetic analysis.

## partitions

You can change the partition which is written, by using the --partition option.
`ppanggolin msa -p pangenome.h5 --partition persistent` for example will compute MSA for all the persistent gene families.

Supported partitions are core, persistent, shell, cloud, softcore, accessory. If you wish to have additional filters, you can raise an issue with your demand, or write a PR directly, most possibilites should be quite straightforward to add.

## source

You can specify whether to use dna or protein sequences for the MSA by using --source. It uses protein sequences by default.

`ppanggolin msa -p pangenome.h5 --source dna`

## phylo

It is also possible to write a single whole genome MSA file, which many phylogenetic softwares accept as input, by using the --phylo option as such:

`ppanggolin msa -p pangenome.h5 --phylo`

This will contatenate all of the family MSA into a single MSA, with one sequence for each genome.

# Info

When computing a pangenome, all of the information about it is saved in the .h5 file, notably parameters used at each step and metrics about the pangenome. You can easily retrieve those informations using the 'info' module.

This command prints information on stdout, and does not write any file.

## Content

This option indicates the following metrics about your pangenome, if they have been computed:
* The total number of genes
* The number of genomes
* The number of gene families
* The number of edges in the pangenome graph
* The number of persistent genes, with the minimal, maximal, sd and mean presence thresholds of the families in this partition
* The number of shell genes, with the minimal, maximal, sd and mean presence thresholds of the families in this partition
* The number of cloud genes, with the minimal, maximal, sd and mean presence thresholds of the families in this partition
* The number of partitions

Additionally, if you have predicted RGPs and spots (with the subcommands 'panrgp', 'rgp' and 'spot', or 'all'), you will have the following metrics:
* The number of RGPs (Regions of Genomic Plasticity)
* The number of spots of insertion

Additionally, if you have predicted modules (with the subcommands 'panmodule', 'module' or 'all'):
* The number of modules
* The number of gene families in modules

It is used as such:
` ppanggolin info -p pangenome.h5 --content` 

## Parameters

This option indicates, for each steps of the analysis, the PPanGGOLiN parameters that were used and the source of the data if appropriate.

It is used as such:

`ppanggolin info -p pangenome.h5 --parameters`

# Metrics
After computing a pangenome, it's interesting to get some metrics about it.
The `metrics` subcommand allow running and compute some analysis and metrics.

All the metrics computed here will be saved in your pangenome file and 
will be easily readable with the `info` subcommand

## Genomic fluidity
The genomic fluidity is described as *a robust metric to categorize the
gene-level similarity among groups of sequenced isolates.* 
[more information here](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-32)

We add the possibility to get genomic fluidity for all the pangenome or 
for specific partition. The genomic fluidity is computable like that :
```
ppanggolin metrics -p pangenome --genome_fluidity
...
Genomes fluidity : all=0.026, shell=0.477, cloud=0.045, accessory=0.554
```
*all* correspond to all the family in the pangenome (core and accessory)

## Module information
It could be necessary to get more information about the modules. 
Here we provide information about families, and we separate modules in 
function of the partition. You can get this supplementary information 
as such :
```
ppanggolin metrics -p pangenome.h5 --info_modules
...
Modules : 3
Families in Modules : 22  (min : 5, max : 9, sd : 2.08, mean : 7.33)
	Sheel specific : 36.36  (sd : 4.62, mean : 2.67)
	Cloud specific : 63.64  (sd : 4.51, mean : 4.67)
```

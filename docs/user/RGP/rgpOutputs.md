## RGP outputs

### RGP

The `plastic_regions.tsv` is a tsv file that lists all of the detected Regions of Genome Plasticity. This requires to have run the RGP detection analysis by either using the `panrgp` command or the `rgp` command.

It can be written with the following command:
`ppanggolin write_pangenome -p pangenome.h5 --regions`

The file has the following format :

| Column | Description |
|--------|-------------|
|region| A unique identifier for the region. This is usually built from the contig it is on, with a number after it.|
|genome| The genome it is in. This is the genome name provided by the user.|
|start| The start position of the RGP in the contig.|
|stop| The stop position of the RGP in the contig.|
|genes| The number of genes included in the RGP.|
|contigBorder| This is a boolean column. If the RGP is on a contig border it will be True, otherwise, it will be False. This often can indicate that, if an RGP is on a contig border it is probably not complete.|
|wholeContig| This is a boolean column. If the RGP is an entire contig, it will be True, and False otherwise. If a RGP is an entire contig it can possibly be a plasmid, a region flanked with repeat sequences or a contaminant.|

### Spots

The `spots.tsv` is a tsv file with two column. It links the spots of `summarize_spots.tsv` with the RGPs of `plastic_regions.tsv`.

It is written with the following command:
`ppanggolin write_pangenome -p pangenome.h5 --spots`

|Column|Description|
|------|------------|
|spot_id| The spot identifier (found in the 'spot' column of `summarize_spots.tsv`).|
|rgp_id| The RGP identifier (found in 'region' column of `plastic_regions.tsv`).|

### Summarize spots

The `summarize_spots.tsv` file is a tsv file that will associate each spot with multiple metrics that can indicate the dynamic of the spot.

It is written with the following command:
`ppanggolin write_pangenome -p pangenome.h5 --spots`

|Column|Description|
|-------|------------|
|spot|The spot identifier. It is unique in the pangenome.|
|nb_rgp|The number of RGPs present in the spot.|
|nb_families| The number of different gene families that are found in the spot.|
|nb_unique_family_sets| The number of RGPs with different gene family content. If two RGPs are identical, they will be counted only once. The difference between this number and the one provided in 'nb_rgp' can be a strong indicator on whether their is a high turnover in gene content in this area or not.|
|mean_nb_genes| The mean number of genes on RGPs in the spot.|
|stdev_nb_genes| The standard deviation of the number of genes in the spot.|
|max_nb_genes| The longest RGP in number of genes of the spot.|
|min_nb_genes| The shortest RGP in number of genes of the spot.|

### Borders

Each spot has at least one set of gene families bordering them. To write the list of gene families bordering a spot, you need to use the following option:
`ppanggolin write_pangenome -p pangenome.h5 --borders`

It will write a .tsv file with 4 columns:

|Column|Description|
|-------|------------|
|spot_id| The spot identifier. It is unique in the pangenome.|
|number| The number of RGPs present in the spot that have those bordering genes.|
|border1| Comma-separated list of gene families of the 1st border.|
|border2| Comma-separated list of gene families of the 2nd border.|

As there can be some variation in the borders, some spots will have multiple borders and as such multiple lines in this file.
The sum of the number for each spot_id should be exactly the number of RGPs in the spot.

## Draw spots

The `draw` command can draw specific spots of interest, whose ID are provided, or all the spots if you wish.
It will also write a gexf file, which corresponds to the gene families and their organization within the spots. It is basically a subgraph of the pangenome, consisting of the spot itself.
The command can be used as such:

`ppanggolin draw -p pangenome.h5 --spots all` will draw an interactive `.html` figure and a `gexf` graph file for all the spots.

If you are interested in only a single spot, you can use its identifier to draw it, as such:

`ppanggolin draw -p pangenome.h5 --spots spot_34` for spot_34, for example.

The interactive figures that are drawn look like this:

```{image} ../../_static/drawspot_example.png
:align: center
```

The plot represents the different gene organizations that are found in the spot. If there are RGPs with identical gene organization, the organization is represented only once (the represented RGP is picked at random among all identical RGPs). The list of RGPs with the same organization is accessible in the file written alongside the figure called `spot_X_identical_rgps.tsv`, with X the spot_id.

They can be edited using the sliders and the radio buttons, to change various graphical parameters, and then the plot itself can be saved using the save button on the right of the screen, if need be.

For the gexf file, you can see how to visualize it in the section about the [pangenome gexf](../PangenomeAnalyses/pangenomeGraphOut.md#pangenome-graph-output).
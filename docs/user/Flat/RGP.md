#### RGP
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

#### Spots

This is a tsv file with two column. It links the spots of 'summarize_spots' with the RGPs of 'plastic_regions'.

It is written with the following command:
`ppanggolin write -p pangenome.h5 --spots`

|column|description|
|------|------------|
|spot_id| The spot identifier (found in the 'spot' column of 'summarize_spots')|
|rgp_id| the RGP identifier (found in 'region' column of 'plastic_regions')|

#### Summarize spots

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

#### Borders

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
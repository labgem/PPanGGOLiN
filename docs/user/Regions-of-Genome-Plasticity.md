From version 1.1.0 and on, it is possible to predict and work with Regions of Genome Plasticity (RGP) using PPanGGOLiN.
RGPs correspond roughly to genomic islands, plasmids, and regions that have been lost in multiple strains. They are areas of the genome where there is a stretch of _shell_ and _cloud_ genes, which can indicate a more plastic area than those made of only _persistent_ genes.

Those analyses can be done using the `ppanggolin panrgp` command directly from your .fasta instead of the `ppanggolin workflow` command. This will write additional files that are related to regions of genome plasticity. Descriptions of those files can be found in the [output wiki](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#plastic-regions).

This analysis can also be run with dedicated subcommands.

# RGP

Once partitions have been computed, you can predict the regions of genome plasticity. 
This subcommand's options are about tuning parameters for the analysis. Details about each parameter can be found in the related [article](https://doi.org/10.1093/bioinformatics/btaa792).

You can do it as such:

`ppanggolin rgp -p pangenome.h5`

This will predict RGPs and store results in the HDF5 file. If you want a list of RGPs for each genome, you can use `ppanggolin write -p pangenome.h5 --regions --output MYOUTPUTDIR`. It will provide the file 'plastic regions' whose format is described [here](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#plastic-regions)

# Spots of insertion

To study RGPs that are found in the same area in different genomes, we gather them into 'spots of insertion'. Those spots are groups of RGPs that do not necessarily have the same gene content but have similar bordering _persistent_ genes. We run those analyses to study the dynamic of gene turnover of large regions in bacterial genomes. In this way, spots of the same pangenome can be compared and their dynamic can be established by comparing their different metrics together. Those metrics are described in the [output wiki](https://github.com/labgem/PPanGGOLiN/wiki/Outputs).

Spots can be computed once RGPs have been predicted. You can do that using:

`ppanggolin spot -p pangenome.h5`

For versions between 1.1.0 and 1.2.12, you can use additional option '--draw_hotspots' which uses [genoplotR](http://genoplotr.r-forge.r-project.org/) to draw those spots in png figures. For versions above 1.2.12, you can use the dedicated subcommand [draw](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#draw), which uses the python library [bokeh](http://docs.bokeh.org/en/latest/) to draw interactive figures which can be visualized and modified directly in the browser.

Information about spots can then be written using `ppanggolin write -p pangenome --spots` which will provide a [file linking RGPs with their spots](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#spots) and a [file showing multiple metrics for each spot](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#summarize-spots)



# RGP cluster based on their gene families

To cluster RGPs (Regions of Genome Plasticity) based on their gene families, you can use the command `panggolin rgp_cluster`.
The panggolin rgp_cluster command performs the following steps to cluster RGPs (Regions of Genome Plasticity) based on their gene families:

1. Calculation of GRR (Gene Repertoire Relatedness): The command calculates the GRR values for all pairs of RGPs. The GRR metric evaluates the similarity between two RGPs by assessing their shared gene families.
2. Graph Construction: The command constructs a graph representation of the RGPs, where each RGP is represented as a node in the graph. The edges between the nodes are weighted using the GRR values, indicating the strength of the relationship between the RGPs.
3. Filtering GRR Values: GRR values below the `--grr_cutoff` threshold (default 0.8) are filtered out to remove noise from the analysis.
4. Louvain Communities Clustering: The Louvain communities clustering algorithm is then applied to the graph. This algorithm identifies clusters of RGPs with similar gene family relationships.

There are three modes available for calculating the GRR value: `min_grr`, `max_grr`, or `incomplete_aware_grr`.
- `min_grr` mode: This mode computes the number of gene families shared between two RGPs and divides it by the smaller number of gene families among the two RGPs.
- `max_grr` mode: In this mode, the number of gene families shared between two RGPs is calculated and divided by the larger number of gene families among the two RGPs.
- `incomplete_aware_grr` (default) mode: If at least one RGP is considered incomplete, which typically happens when it is located at the border of a contig, the `min_grr` mode is used. Otherwise, the `max_grr` mode is applied. This mode is useful to correctly cluster incomplete RGPs.


The resulting RGP clusters are stored in a tsv file with the folowing columns:

| column  | description                  |
|---------|------------------------------|
| RGP     | The unique region identifier |
| cluster | The cluster id of the RGP    |
| spot_id    | the spot ID of the RGP       |



The command also generates an RGP graph in the gexf format, which can be utilized to explore the RGP clusters along with their spots of insertion. In this graph identical RGPs with the same family content and with the same spot are merged into a single node to simplify the graph representation. This feature can be disable with the parameter `--no_identical_rgp_merging`.


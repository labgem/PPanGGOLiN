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

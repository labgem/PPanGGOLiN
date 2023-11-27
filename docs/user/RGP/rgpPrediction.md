## PanRGP


This command works exactly like 'workflow'. The difference is that it will run more analysis related to [Regions of Genome Plasticity](#RGP-section).
You can use the panrgp command as follow:

```bash
ppanggolin panrgp --fasta ORGANISMS_FASTA_LIST
```

The rgp analysis is launched after the pangenome partitionning and use the default parameters. 
If you want to tune the rgp detection, you can use the `rgp` command after the `workflow` command.


More detail about RGP detection [here](#RGP-section) and in the [panRGP publication](https://doi.org/10.1093/bioinformatics/btaa792)


## RGP detection

Once partitions have been computed, you can predict the regions of genome plasticity. 
This subcommand's options are about tuning parameters for the analysis. Details about each parameter can be found in the related [article](https://doi.org/10.1093/bioinformatics/btaa792).

You can do it as such:

`ppanggolin rgp -p pangenome.h5`

This will predict RGPs and store results in the HDF5 file. If you want a list of RGPs for each genome, you can use `ppanggolin write -p pangenome.h5 --regions --output MYOUTPUTDIR`. It will provide the file 'plastic regions' whose format is described [here](rgpOutputs.md#RGP)

## Spot prediction


To study RGPs that are found in the same area in different genomes, we gather them into 'spots of insertion'. Those spots are groups of RGPs that do not necessarily have the same gene content but have similar bordering _persistent_ genes. We run those analyses to study the dynamic of gene turnover of large regions in bacterial genomes. In this way, spots of the same pangenome can be compared and their dynamic can be established by comparing their different metrics together. Those metrics are described in the [RGP and spot output section]g(rgpOutputs.md#Spots).

Spots can be computed once RGPs have been predicted. You can do that using:

`ppanggolin spot -p pangenome.h5`

For versions between 1.1.0 and 1.2.12, you can use additional option `--draw_hotspots` which uses [genoplotR](http://genoplotr.r-forge.r-project.org/) to draw those spots in png figures. For versions above 1.2.12, you can use the dedicated subcommand `draw` with the argument `--spots`, to draw interactive figures with the python library [bokeh](http://docs.bokeh.org/en/latest/)  which can be visualized and modified directly in the browser. This plot is discribe [here](rpOutputs.md#draw-spots)

Information about spots can then be written using `ppanggolin write -p pangenome --spots` which will provide a [file linking RGPs with their spots](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#spots) and a [file showing multiple metrics for each spot](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#summarize-spots)
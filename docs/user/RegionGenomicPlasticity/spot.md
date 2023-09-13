
To study RGPs that are found in the same area in different genomes, we gather them into 'spots of insertion'. Those spots are groups of RGPs that do not necessarily have the same gene content but have similar bordering _persistent_ genes. We run those analyses to study the dynamic of gene turnover of large regions in bacterial genomes. In this way, spots of the same pangenome can be compared and their dynamic can be established by comparing their different metrics together. Those metrics are described in the [output wiki](https://github.com/labgem/PPanGGOLiN/wiki/Outputs).

Spots can be computed once RGPs have been predicted. You can do that using:

`ppanggolin spot -p pangenome.h5`

For versions between 1.1.0 and 1.2.12, you can use additional option '--draw_hotspots' which uses [genoplotR](http://genoplotr.r-forge.r-project.org/) to draw those spots in png figures. For versions above 1.2.12, you can use the dedicated subcommand [draw](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#draw), which uses the python library [bokeh](http://docs.bokeh.org/en/latest/) to draw interactive figures which can be visualized and modified directly in the browser.

Information about spots can then be written using `ppanggolin write -p pangenome --spots` which will provide a [file linking RGPs with their spots](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#spots) and a [file showing multiple metrics for each spot](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#summarize-spots)
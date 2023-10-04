
A tile plot is a heatmap representing the gene families (y axis) in the organisms (x axis) making up your pangenome. The tiles on the graph will be colored if the gene family is present in an organism and uncolored if absent. The gene families are ordered by partition, and the genomes are ordered by a hierarchical clustering based on their shared gene families (basically two genomes that are close together in terms of gene family composition will be close together on the figure).

This plot is quite helpful to observe potential structures in your pangenome, and can also help you to identify eventual outliers. You can interact with it, and mousing over a tile in the plot will indicate to you which is the gene identifier(s), the gene family and the organism that corresponds to the tile.

If you build your pangenome using the 'workflow' subcommand and you have more than 500 organisms, only the 'shell' and the 'persistent' partitions will be drawn, leaving out the 'cloud' as the figure tends to be too heavy for a browser to open it otherwise.

It can be generated using the 'draw' subcommand as such : 

`ppanggolin draw -p pangenome.h5 --tile_plot`


```{image} ../_static/tile_plot.png
:align: center
```

and if you do not want the 'cloud' gene families as it is a lot of data and can be hard to open with a browser sometimes, you can use the following option : 

`ppanggolin draw -p pangenome.h5 --tile_plot --nocloud`

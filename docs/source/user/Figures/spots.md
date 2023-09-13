
For versions 1.2.30 and above, the 'draw' command can draw specific spots of interest, whose ID are provided, or all the spots if you wish.
It will also write a gexf file, which corresponds to the gene families and their organization within the spots. It is basically a subgraph of the pangenome, consisting of the spot itself.
The command can be used as such:

`ppanggolin draw -p pangenome.h5 --spots all` will draw an interactive .html figure and a gexf file for all the spots.

If you are interested in only a single spot, you can use its identifier to draw it, as such:

`ppanggolin draw -p pangenome.h5 --spots spot_34` for spot_34, for example.

The interactive figures that are drawn look like this:

```{image} ../_static/drawspot_example.png
:align: center
```

The plot represents the different gene organizations that are found in the spot. If there are RGPs with identical gene organization, the organization is represented only once (the represented RGP is picked at random among all identical RGPs). The list of RGPs with the same organization is accessible in the file written alongside the figure called 'spot_X_identical_rgps.tsv', with X the spot_id.

They can be edited using the sliders and the radio buttons, to change various graphical parameters, and then the plot itself can be saved using the save button on the right of the screen, if need be.
For the gexf file, you can see how to visualize it in the section about the [pangenome gexf](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#gexf-and-light-gexf).
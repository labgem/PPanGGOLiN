(output)=
# PPanGGOLiN outputs 

PPanGGOLiN provides multiple outputs to describe a pangenome. In this section the different outputs will be described.

In most cases it will provide with a HDF-5 file named "pangenome.h5". This file stores all the information about your pangenome and the analysis that were run. If given to ppanggolin through most of the subcommands, it will read information from it. This is practical as you can regenerate figures or output files, or rerun parts of the analysis without redoing everything.

In this section, each parts will describe a possible output of PPanGGOLiN, and will be commented with the command line that generates it using the HDF5 file, which is assumed to be called 'pangenome.h5'.

When using the same subcommand (like 'write' or 'draw' that can help you generate multiple file each), you can provide multiple options to write all of the file formats that you desire at once.

## PPanGGOLiN figures outputs

### U-shaped plot
```{include} Figures/Uplot.md
```

### tile plot
```{include} Figures/tilePlot.md
```

### Spot plots
```{include} Figures/spots.md
```

## Rarefaction
```{include} Figures/rarefaction.md
```

## Write
### Organisms statistics
```{include} Flat/orgStat.md
```

### pangenomeGraph files
The pangenome's graph can be given through multiple data formats, in order to manipulate it with other softwares.

#### gexf and light gexf
```{include} graphOut/GEXF.md
```


#### json
```{include} graphOut/JSON.md
```

```{include} Flat/presAbs.md
```

### mean persistent duplication
```{include} Flat/dupplication.md
```

### partitions
```{include} Flat/partition.md
```

### projection
```{include} Flat/projection.md
```

### Gene families and genes
```{include} Flat/fam2gen.md
```

### Genomic Island
```{include} Flat/RGP.md
```

### Modules
```{include} Flat/module.md
```

## Fasta
```{include} sequence/fasta.md
```

## MSA
```{include} sequence/MSA.md
```

## Info
```{include} Flat/info.md
```

## Metrics
```{include} Flat/metrics.md
```
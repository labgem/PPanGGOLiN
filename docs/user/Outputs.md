(output)=
# PPanGGOLiN outputs 

PPanGGOLiN provides multiple outputs to describe a pangenome. In this section the different outputs will be described.

In most cases it will provide with a HDF-5 file named "pangenome.h5". This file stores all the information about your pangenome and the analysis that were run. If given to ppanggolin through most of the subcommands, it will read information from it. This is practical as you can regenerate figures or output files, or rerun parts of the analysis without redoing everything.

In this section, each part will describe a possible output of PPanGGOLiN, and will be commented with the command line that generates it using the HDF5 file, which is assumed to be called 'pangenome.h5'.

When using the same subcommand (like 'write_pangenome' or 'draw' that can help you generate multiple file each), you can provide multiple options to write all of the file formats that you desire at once.

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

### Rarefaction
```{include} Figures/rarefaction.md
```

##  Write flat outputs describing the pangenome

Writes 'flat' files that describe the pangenome and its elements with the command `write_pangenome`.

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

### Gene families and genes
```{include} Flat/fam2gen.md
```

### Genomic Island
```{include} Flat/RGP.md
```

### Modules
```{include} Flat/module.md
```

### Partitions
```{include} Flat/partition.md
```

## Write genome files with pangenome annotations
The `write_genomes` command creates 'flat' files representing genomes with their pangenome annotations.

To generate output exclusively for particular genomes, users can utilize the `--organisms` argument. This argument accepts a list of organism names, either directly entered in the command line (comma-separated) or referenced from a file where each line contains a single organism name.


### Genes table with pangenome annotations
```{include} Flat/tables.md
```
### GFF file
```{include} Flat/gff.md
```
### JSON Map for Proksee visualisation
```{include} Flat/proksee.md
```
### Adding Fasta Sequences into GFF and proksee JSON map Files

```{include} Flat/genomes_fasta.md
```

### Incorporating Metadata into Tables, GFF, and Proksee Files
```{include} Flat/genomes_metadata.md
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
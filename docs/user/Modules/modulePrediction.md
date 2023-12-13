<!-- # Conserved module prediction -->

PPanGGOLiN can predict and work with conserved modules, which are groups of genes that are part of the variable genome, and often found together across the genomes of the pangenome. These conserved modules may also be potential functional modules.

Further details can be found in the [panModule preprint](https://doi.org/10.1101/2021.12.06.471380)

## The panModule workflow

The panModule workflow facilitates the generation of a pangenome with predicted conserved modules from a specified set of genomes. This command extends the functionality of the `workflow` command by detecting conserved modules. Additionally, it generates descriptive tsv files detailing the predicted modules, whose format are detailed [here](./moduleOutputs.md).

To execute the panModule workflow, use the following command: 

```bash
ppanggolin panmodule --fasta GENOME_LIST_FILE
```
Replace `GENOME_LIST_FILE` with a tab-separated file listing the genome names, and the fasta file path of their genomic sequences as described [here](../PangenomeAnalyses/pangenomeAnnotation.md#annotate-from-fasta-files). Alternatively, you can provide a list of GFF/GBFF files as input by using the `--anno` parameter, similar to how it is used in the workflow and annotate commands.

The panmodule workflow predicts modules using default parameters. To fine-tune the detection, you can use the `module` command on a partioned pangenome acquired through the `workflow` for example or use a configuration file, as described [here](../practicalInformation.md#configuration-file). 


## Predict conserved module

The `module` command predicts conserved modules on an partioned pangenome. The command has several options for tuning the prediction. Details about each parameter are available in the related [preprint](https://www.biorxiv.org/content/10.1101/2021.12.06.471380v1).

The command can be used simply as such:

```bash
ppanggolin module -p pangenome.h5
```

This will predict modules and store the results in the HDF5 pangenome file. If you wish to have descriptive tsv files, whose format is detailed [here](./moduleOutputs.md), you can use the `write_pangenome` command with the flag `--modules`:
```bash
ppanggolin write_pangenome -p pangenome.h5 --modules --output MYOUTPUTDIR
```

If spots of insertion have been predicted in you pangenome using the `spot` command (or inside the `panrgp` or `all` workflow commands), you can also list the associations between the predicted spots and the predicted modules as such:

```bash
ppanggolin write_pangenome -p pangenome.h5 --spot_modules --output MYOUTPUTDIR
```


The format of each file is given [here](./moduleOutputs.md)

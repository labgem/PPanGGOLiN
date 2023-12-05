# Conserved module prediction

PPanGGOLiN is able to predict and work with conserved modules. Modules are groups of genes that are part of the variable genome, and often found together across the genomes of the pangenome. As such, they are conserved modules and potential functional modules.

Further details can be found in the [panModule preprint](https://doi.org/10.1101/2021.12.06.471380)

## The panModule workflow

The panModule workflow facilitates the generation of a pangenome with predicted conserved modules from a specified set of genomes. This command extends the functionality of the `workflow` command by detecting conserved modules. Additionally, it generates descriptive TSV files detailing the predicted modules, whose format are detailed [here](./moduleOutputs.md).

To execute the panModule workflow, use the following command: 

```bash
ppanggolin panmodule --fasta GENOME_LIST_FILE
```
Replace `GENOME_LIST_FILE` with a tab-separated file listing the genome names, and the fasta filepath of its genomic sequences as described [here](../PangenomeAnalyses/pangenomeAnnotation.md#annotate-fasta-file). Alternatively, provide a list of GFF/GBFF files as input by utilizing the `--anno` parameter, similar to how it's used in the workflow and annotate commands.

The panmodule workflow predict modules with default parameters. To fine-tune the detection, you have the option to use the `module` command on a partionned pangenome acquired through the `workflow` for example or use a configuration file, detailed further [here](../practicalInformation.md#configuration-file). 


## Predict conserved module

The `module` command predicts conserved modules on an partioned pangenome. The command has several options for tuning the prediction. Details about each parameter are available in the related [preprint](https://www.biorxiv.org/content/10.1101/2021.12.06.471380v1).

The command can be used simply as such:

`ppanggolin module -p pangenome.h5`

This will predict modules and store the results in the HDF5 pangenome file. If you wish to have descriptive tsv files, whose format is detailed [here](./moduleOutputs.md), you can use the `write_pangenome` command with the flag `--modules`:

`ppanggolin write_pangenome -p pangenome.h5 --modules --output MYOUTPUTDIR`.

If spots of insertion have been predicited in you pangenome using the `spot` command (or inside the `panrgp` or `all` workflow commands), you can also list the associations between the predicted spots and the predicted modules as such:

`ppanggolin write_pangenome -p pangenome.h5 --spot_modules --output MYOUTPUTDIR`

The format of each file is given [here](./moduleOutputs.md)
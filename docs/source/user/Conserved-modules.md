(module-section)=
# Conserved Module
From version 1.2.0, it is possible to predict and work with conserved modules using PPanGGOLiN. Modules are groups of genes that are part of the variable genome, and often found together in the different genomes. As such, they are conserved modules and potential functional modules.

Once partitions have been computed, you can predict conserved modules. All the options of the `module` subcommand are for tuning the parameters for the analysis.
Details about each parameter and what they do is available in the related [preprint](https://www.biorxiv.org/content/10.1101/2021.12.06.471380v1).

The command can be used simply as such:

`ppanggolin module -p pangenome.h5`

This will predict modules and store the results in the HDF5 file. If you wish to have descriptive tsv files, whose format is detailed [here](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#functional-modules), you can use:

`ppanggolin write -p pangenome.h5 --modules --output MYOUTPUTDIR`.

If your pangenome has spots of insertion that were predicted using the `spot` command (or the `panrgp` or `all` commands), you can also list the associations between the predicted spots and the predicted modules as such:

`ppanggolin write -p pangenome.h5 --spot_modules --output MYOUTPUTDIR`

The format of each file is given [here](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#spot-modules)
### gene presence absence

This file is basically a presence absence matrix. The columns are the genomes used to build the pangenome, the lines are the gene families. The identifier of the gene family is the gene identifier chosen as a representative.
 There is a 1 if the gene family is present in a genome, and 0 otherwise. It follows the exact same format than the 'gene_presence_absence.Rtab' file that you get from the pangenome analysis software [roary](https://sanger-pathogens.github.io/Roary/)

It can be generated using the 'write' subcommand as such : 

`ppanggolin write_pangenome -p pangenome.h5 --Rtab`

### matrix

This file is a .csv file following a format alike the gene_presence_absence.csv file generated by [roary](https://sanger-pathogens.github.io/Roary/), and works with [scoary](https://github.com/AdmiralenOla/Scoary) if you want to do pangenome-wide association studies.

It can be generated using the 'write' subcommand as such : 

`ppanggolin write_pangenome -p pangenome.h5 --csv`
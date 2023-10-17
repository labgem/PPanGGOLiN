
You can write a list containing the gene family assigned to every single gene of your pangenome, in a file format extactly like the one provided by [MMseqs2](https://github.com/soedinglab/MMseqs2) through its subcommand 'createtsv'.
It is basically a three-column file listing the gene family name in the first column, and the gene names in the second. A third column is either empty, or has an "F" in it. In that case it indicates that the gene is potentially a gene fragment and not complete. This will be indicated only if the [defragmentation](https://github.com/labgem/PPanGGOLiN/wiki/PPanGGOLiN---step-by-step-pangenome-analysis#defragmentation) pipeline is used.

You can obtain it as such :  

`ppanggolin write -p pangenome.h5 --families_tsv`
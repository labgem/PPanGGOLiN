This option writes in a 'tables' directory. There will be a file written in the .tsv file format for every single genome in the pangenome.
The columns of this file are described in the following table : 

| Column               | Description                                                                                                                    |
|----------------------|--------------------------------------------------------------------------------------------------------------------------------|
| gene                 | the unique identifier of the gene                                                                                              |
| contig               | the contig that the gene is on                                                                                                 |
| start                | the start position of the gene                                                                                                 |
| stop                 | the stop position of the gene                                                                                                  |
| strand               | The strand that the gene is on                                                                                                 |
| ori                  | Will be T if the gene name is dnaA                                                                                                              |
| family               | the family identifier to which the gene belongs to                                                                             |
| nb_copy_in_org       | The number of copy of the family in the organism (basically, if 1, the gene has no closely related paralog in that organism) |
| partition            | the partition to which the gene family of the gene belongs to                                                                  |
| persistent_neighbors | The number of neighbors classified as 'persistent' in the pangenome graph                                                      |
| shell_neighbors      | The number of neighbors classified as 'shell' in the pangenome graph                                                           |
| cloud_neighbors      | The number of neighbors classidied as 'cloud' in the pangenome graph                                                           |

Those files can be generated as such : 

`ppanggolin write_genomes -p pangenome.h5 --tables`
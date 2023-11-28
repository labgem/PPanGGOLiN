Those files will be stored in the 'partitions' directory and will be named after the partition that they represent (like persistent.txt for the persistent partition). In each of those file there will be a list of gene family identifiers that correspond to the gene families belonging to that partition, one family per line, should you need it for your pipelines or during your analysis.

You can generate those files as such :  

` ppanggolin write_pangenome -p pangenome.h5 --partitions`
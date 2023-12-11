
This step will assign gene families either to the 'persistent', 'shell', or 'cloud' partitions. 

The 'persistent' partition includes gene families with genes commonly found throughout the species
It corresponds to essential genes, genes required for important metabolic pathways and genes that define the metabolic and biosynthetic capabilities of the taxonomic group.

The 'shell' partition includes gene families with genes found in some individuals. 
These genes, frequently acquired via horizontal gene transfers, typically encode functions related to environmental adaptation, pathogenicity, virulence, or the production of secondary metabolites

The 'cloud' partition includes gene families with rare genes found in one, or very few, individuals. 
Most of the genes were associated with phage-related genes. 
They probably all were acquired through horizontal gene transfers. 
Antibiotic resistance genes were often found to be belonging to the cloud genome, as well as plasmid genes.

It can be realized through the following subcommand : 

`ppanggolin partition -p pangenome.h5`

The command also has quite a few options. 
Most of them are not self-explanatory. 
If you want to know what they do, you should read the PPanGGOLiN paper (you can read it [here](https://journals.plos.org/ploscompbiol/article?rev=2&id=10.1371/journal.pcbi.1007732)) where the statistical methods used are thoroughly described.

The one parameter that might be of importance is the '-K', or '--nb_of_partitions' parameter. 
This will define the number of classes (`K`) used to partition the pangenome. 
If you anticipate well-defined subpopulations within your pangenome and know their exact number, this approach can be particularly useful. For metagenome-assembled genomes (MAGs), which are often incomplete, it is typically advised to set a fixed value of K=3.
If the number of subpopulations is unknown, it will be automatically determined using the ICL criterion.
The idea is that the most present partition will be 'persistent', the least present will be 'cloud', and all the others will be 'shell'. 
The number of partitions corresponding to the shell will be the number of expected subpopulations in your pangenome. 
(for instance, if you expect 5 subpopulations, you could use -K 7). 

In most cases, you should let the statistical criterion used by PPanGGOLiN find the optimal number of partitions for you.

All the results will be added to the given 'pangenome.h5' input file.

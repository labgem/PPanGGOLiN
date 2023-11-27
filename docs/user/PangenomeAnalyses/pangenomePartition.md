
This is the step that will assign gene families to the 'persistent', 'shell', or 'cloud' partitions. 


The 'persistent' partition will group genes that are present throughout the entire species. 
They will be essential genes, genes required for important metabolic pathways and genes that define the metabolic and biosynthetic capabilities of the taxonomic group.

The 'shell' partition groups genes that are present in only some individuals. 
Those are often genes that were acquired through horizontal gene transfers and encode for functions involved in environmental adaptations, pathogenicity, virulence or encoding secondary metabolites for example.

The 'cloud' partition groups genes that are very rare in the pangenome and found in one, or very few, individuals. 
Most of the genes were associated with phage-related genes. 
They probably all were acquired through horizontal gene transfers. 
Antibiotic resistance genes were often found to be belonging to the cloud genome, as well as plasmid genes.

It can be realized through the following subcommand : 

`ppanggolin partition -p pangenome.h5`

It also has quite a few options. 
Most of them are not self-explanatory. 
If you want to know what they do, you should read the PPanGGOLiN paper (you can read it [here](https://journals.plos.org/ploscompbiol/article?rev=2&id=10.1371/journal.pcbi.1007732)) where the statistical methods used are thoroughly described.

The one parameter that might be of importance is the '-K', or '--nb_of_partitions' parameter. 
This will define the number of classes used to partition the pangenome. 
This may be of use if you expect to have well-defined subpopulations in your pangenome, and you know exactly how many. 
If not, that number is detected automatically through an ICL criterion. 
The idea is that the most present partition will be 'persistent', the least present will be 'cloud', and all the others will be 'shell'. 
The number of partitions corresponding to the shell will be the number of expected subpopulations in your pangenome. 
(So if you expect 5 subpopulations, you could use -K 7). 


In most cases, you should let the statistical criterion used by PPanGGOLiN find the optimal number of partitions for you.

All the results will be added to the given 'pangenome.h5' input file.
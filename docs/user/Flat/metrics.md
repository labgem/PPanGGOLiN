After computing a pangenome, it's interesting to get some metrics about it.
The `metrics` subcommand allow running and compute some analysis and metrics.

All the metrics computed here will be saved in your pangenome file and 
will be easily readable with the `info` subcommand

### Genomic fluidity
The genomic fluidity is described as *a robust metric to categorize the
gene-level similarity among groups of sequenced isolates.* 
[more information here](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-32)

We add the possibility to get genomic fluidity for all the pangenome or 
for specific partition. The genomic fluidity is computable like that :
```
ppanggolin metrics -p pangenome --genome_fluidity
...
Genomes fluidity : all=0.026, shell=0.477, cloud=0.045, accessory=0.554
```
*all* correspond to all the family in the pangenome (core and accessory)

### Module information
It could be necessary to get more information about the modules. 
Here we provide information about families, and we separate modules in 
function of the partition. You can get this supplementary information 
as such :
```
ppanggolin metrics -p pangenome.h5 --info_modules
...
Modules : 3
Families in Modules : 22  (min : 5, max : 9, sd : 2.08, mean : 7.33)
	Sheel specific : 36.36  (sd : 4.62, mean : 2.67)
	Cloud specific : 63.64  (sd : 4.51, mean : 4.67)
```

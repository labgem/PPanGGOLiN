### Pangenome metrics


#### Genomic fluidity

The genomic fluidity is described as *a robust metric to categorize the
gene-level similarity among groups of sequenced isolates.* 
[more information here](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-32)

We add the possibility to get genomic fluidity for all the pangenome or 
for specific partition. The genomic fluidity is computable like that (with below an example result):
```bash
ppanggolin metrics -p pangenome --genome_fluidity
...
Genomes fluidity: all=0.026, shell=0.477, cloud=0.045, accessory=0.554
```
*all* correspond to all the family in the pangenome (core and accessory)
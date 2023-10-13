This command is available from 1.1.98 and on.
This command can be used to write fasta sequences of the pangenome or of parts of the pangenome. Most options expect a partition to write. Available partitions are:
* 'all' for the entire pangenome.
* 'Persistent' for persistent families
* 'Shell' for shell genes or families
* 'Cloud' for cloud genes or families
* 'rgp' for genes or families found in RGPs
* 'core' for core genes or families
* 'softcore' for softcore genes or families

When using the 'softcore' filter, the '--soft_core' option can be used to modily the threshold used to determine what is part of the soft core. It is at 0.95 by default.

### Genes

This option can be used to write the nucleotide CDS sequences. It can be used as such, to write all of the genes of the pangenome for example:

```ppanggolin fasta -p pangenome.h5 --output MY_GENES --genes all```

Or to write only the persistent genes:

```ppanggolin fasta -p pangenome.h5 --output MY_GENES --genes persistent```


### Protein families

This option can be used to write the protein sequences of the representative sequences for each family. It can be used as such for all families:

```ppanggolin fasta -p pangenome.h5 --output MY_PROT --prot_families all```

or for all of the shell families for example:

```ppanggolin fasta -p pangenome.h5 --output MY_PROT --prot_families shell```


### Gene families

This option can be used to write the gene sequences of the representative sequences for each family. It can be used as such:

```ppanggolin fasta -p pangenome.h5 --output MY_GENES --gene_families all```

or for the cloud families for example:

```ppanggolin fasta -p pangenome.h5 --output MY_GENES --gene_families cloud```

### Regions

This option can be used to write the nucleotide sequences of the detected RGPs. It requires the fasta sequences that were used to compute the pangenome as they were provided originally when you computed your pangenome. This command only has two filters:
* all, for all regions
* complete, for only the 'complete' regions which are not on a contig border

It can be used as such:

```ppanggolin fasta -p pangenome.h5 --output MYREGIONS --regions all --fasta organisms.fasta.list```

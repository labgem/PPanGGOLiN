
# Write pangenome sequences

The `fasta` command can be used to write sequences of the pangenome or specific parts of the pangenome in FASTA format. 

Most options require a partition.

Available partitions are:
* `all` for the entire pangenome.
* `Persistent` for persistent families
* `Shell` for shell genes or families
* `Cloud` for cloud genes or families
* `rgp` for genes or families found in RGPs
* `core` for core genes or families
* `softcore` for softcore genes or families

When using the `softcore` filter, the `--soft_core` option can be used to modify the threshold used to determine what is part of the softcore. It is set to 0.95 by default.

## Genes

This option can be used to write the nucleotide CDS sequences. It can be used as such, to write all of the genes of the pangenome for example:

```bash
ppanggolin fasta -p pangenome.h5 --output MY_GENES --genes all
```

Or to write only the persistent genes:

```bash
ppanggolin fasta -p pangenome.h5 --output MY_GENES --genes persistent
```


## Protein families

This option can be used to write the protein sequences of the representative sequences for each family. It can be used as such for all families:

```bash
ppanggolin fasta -p pangenome.h5 --output MY_PROT --prot_families all
```

or for all of the shell families for example:

```bash
ppanggolin fasta -p pangenome.h5 --output MY_PROT --prot_families shell
```


## Gene families

This option can be used to write the gene sequences of the representative sequences for each family. It can be used as such:

```bash
ppanggolin fasta -p pangenome.h5 --output MY_GENES_FAMILIES --gene_families all
```

or for the cloud families for example:

```bash
ppanggolin fasta -p pangenome.h5 --output MY_GENES_FAMILIES --gene_families cloud
```

## Regions

This option can be used to write the nucleotide sequences of the detected RGPs.
It requires the fasta sequences used to compute the pangenome, as originally provided when you computed your pangenome.

This command has only two filters:
* all, for all regions
* complete, for only the 'complete' regions which are not on a contig border

It can be used as such:

```bash
ppanggolin fasta -p pangenome.h5 --output MY_REGIONS --regions all --fasta genomes.fasta.list
```

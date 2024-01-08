<!-- # Metrics -->

After computing a pangenome, it's interesting to get some metrics about it.
The aim of the `metrics` subcommand is to compute comprehensive metrics describing the pangenome.

All the metrics computed are saved in the pangenome file and 
are displayed by the `info` subcommand with the flag `--content`.


### Genomic fluidity

The genomic fluidity is described as *a robust metric to categorize the
gene-level similarity among groups of sequenced isolates.* 
[more information here](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-32)

We have added the possibility to get genomic fluidity for the whole pangenome or 
for a specific partition. Genomic fluidity is computable as follows:

```bash
ppanggolin metrics -p pangenome --genome_fluidity

```

```yaml
Genomes_fluidity:
    all: 0.139
    shell: 0.527
    cloud: 0.385
    accessory: 0.62
```
*all* correspond to all the family in the pangenome (core and accessory)


```{note}
Currently, the `metrics` command only computes fluidity. However, additional metrics may be added in the future. If you have any ideas for metrics that describe the pangenome, please open an issue! 
```


(RGP-section)=
# Regions of Genome Plasticity
From version 1.1.0 and on, it is possible to predict and work with Regions of Genome Plasticity (RGP) using PPanGGOLiN.
RGPs correspond roughly to genomic islands, plasmids, and regions that have been lost in multiple strains. They are areas of the genome where there is a stretch of _shell_ and _cloud_ genes, which can indicate a more plastic area than those made of only _persistent_ genes.

Those analyses can be done using the `ppanggolin panrgp` command directly from your .fasta instead of the `ppanggolin workflow` command. This will write additional files that are related to regions of genome plasticity. Descriptions of those files can be found in the [output wiki](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#plastic-regions).

This analysis can also be run with dedicated subcommands.

## RGP prediction
```{include} RegionGenomicPlasticity/RGP.md
```

## Spots of insertion detection
```{include} RegionGenomicPlasticity/spot.md
```

## RGP clustering based on their gene families
```{include} RegionGenomicPlasticity/RGPclust.md
```
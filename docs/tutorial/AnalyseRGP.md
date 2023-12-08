# Region of genomic plasticity prediction
## Predict the region of genomic plasticity from the partitioned pangenome graph

"Regions of genome plasticity (**RGPs**) are clusters of genes located in highly variable genomic regions. Most of them arise from HGT and correspond to genomic islands (GIs)." [Bazin et al. 2020](https://doi.org/10.1093/bioinformatics/btaa792)

We are going to start from the partitioned pangenome graph build in the [previous step](#build-and-partition-a-pangenome-graph). 

You can predict the RGP with PPanGGOLiN with this command:

```
ppanggolin rgp -p B_janonicum_results/pangenome.h5
```

As explain before, PPanGGOLiN will predict and store the RGPs in the pangenome file. Follow the 

## Analyse RGPs predicted by PPanGGOLiN
### Get the list of RGPs predicted by PPanGGOLiN
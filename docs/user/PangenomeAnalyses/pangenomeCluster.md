### Cluster genes into gene families

Once we have annotated genomes, we need to compare them to know which are similar, and to build gene families through this information. 

If you provided .fasta files or annotation files with gene sequences in them, clustering can be run directly by providing the .h5 file that was generated, as such : 

```
ppanggolin cluster -p pangenome.h5
```

PPanGGOLiN will call [MMseqs2](https://github.com/soedinglab/MMseqs2) to run the clustering on all the protein sequences by searching for connected components for the clustering step. 
You can tune its parameters using `--identity`(default 0.8) and `--coverage`(default 0.8). 
You can use other clustering algorithms of MMseqs by using --mode (default 1). 
Both protein sequences have to be covered by at least the proportion indicated by --coverage.

#### How to customize MMSeqs2 clustering
```{attention}
All the MMSeqs2 options are not available in PPanGGOLiN if you want a complete view of MMSeqs2 option take a look at their documentation and you can provide your custom clustering as describe in the [next part](#read-clustering)
```

[//]: # (TODO complete this part)

(read-clustering)=
### Providing your gene families

If you do not want to use MMseqs2 and provide your clusters (or gene families) you can do so only if you provided the annotations in the first step. 
In the case of gff3 files, the 'ID' field in the 9th column is expected as a gene id. 
In the case of gbff or gbk files, the 'locus_tag' is used as a gene id, except with files coming from MaGe or from SEED, where the id provided in the 'db_xref' field is used.

You will need to provide a .tsv file. 
The first column indicates the cluster id, and the second column indicates a unique gene id that is used in the annotation files. 
There is a single gene id per line.

You can do that through the command line : 

`ppanggolin cluster -p pangenome.h5 --clusters MY_CLUSTERS_FILE`

An example of what MY_CLUSTERS_FILE should look like is provided [here](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/clusters.tsv)


### Defragmentation

We noticed that most of the cloud genes in the pangenome are fragments of 'shell' or 'persistent' genes, and so not informative on the pangenome's diversity. 
We added another workflow to reduce the number of gene families and reduce the computational load by trying to associate fragments to their original gene families.
It adds a step to the clustering described previously. 
It will compare all the gene families representative protein sequences using the same identity threshold as the first step. 
It will also use the same coverage threshold, but only the smallest of both protein sequences have to be covered by at least the value indicated by `--coverage`.

After that, we build a similarity graph where the edges are the hits given by the comparison, and the nodes are the original gene families. 
Then we iterate on all nodes and compare them to their neighbors. 
If the neighbor of a node is more numerous (has more members in the cluster it represents) and its representative sequence is longer, that node (and all the genes associated) is associated with the neighbor. 
The genes associated with this node are defined as 'fragments' of the gene family represented by the longer and more numerous neighboring node.

To avoid using it, you can run the following:

```
ppanggolin cluster -p pangenome.h5 --no_defrag
```

In any case and whichever pipeline you use, in the end, the gene families will be saved in the 'pangenome.h5' given as input.

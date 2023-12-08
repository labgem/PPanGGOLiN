### Cluster genes into gene families
 
After annotating genomes, we must compare them to determine any similarities and build gene families using this information.

If .fasta files or annotation files containing gene sequences were provided, clustering can be executed by using the generated .h5 file, as such:

```
ppanggolin cluster -p pangenome.h5
```

#### How to customize MMSeqs2 clustering
```{warning}
Not all MMSeqs2 options are available in PPanGGOLiN. For a comprehensive overview of MMSeqs2 options, please refer to their documentation. To create your own custom clustering, please follow the instructions detailed in the [dedicated section](#read-clustering).
```

PPanGGOLiN will run [MMseqs2](https://github.com/soedinglab/MMseqs2) to perform clustering on all the protein sequences by searching for connected components for the clustering step.

##### How to set the identity and coverage parameters

PPanGGOLiN enables the setting of two essential parameters for gene clustering: **identity** and **coverage**. These parameters can be easily adjusted using `--identity` (default: 0.8) and `--coverage` (default: 0.8). The default values were selected as they are the most widely used and effective parameters for aligning and clustering sequences at the species level.
 
Be aware that if you decrease identity and/or coverage, more genes will be clustered together in the same family. 
This will ultimately decrease the number of families and affect all subsequent steps.  

```{note}
The chosen coverage mode in PPanGGOLiN requires both protein sequences to be covered by at least the proportion specified by --coverage.
```

##### How to set the clustering mode

MMSeqs provides 3 different [clustering mode](https://github.com/soedinglab/MMseqs2/wiki#clustering-modes).
By default the clustering mode utilises the _single linkage_ (or _connected component_) algorithm.

Another option is the _set cover_ algorithm, which can be employed using `--mode 1`.

Additionally, the clustering algorithms of MMseqs, similar to CD-Hit, can be selected with `--mode 2` or its low memory version through `--mode 3`.

(read-clustering)=
### Providing your gene families
 
If you want to provide your own clusters (or gene families), you must have provided the annotations in the first step. 
For gff3 files, the expected gene id is the 'ID' field in the 9th column. 
In the case of gbff or gbk files, use 'locus_tag' as a gene id, unless you are working with files from MaGe/MicroScope or SEED, where the id in the 'db_xref' field is used instead.

You will need to provide a .tsv file with a single gene id per line.
The first column should indicate the cluster id, and the second column should indicate the unique gene id used in the annotation files.

You can do this from the command line: 

`ppanggolin cluster -p pangenome.h5 --clusters MY_CLUSTERS_FILE`

An example of what MY_CLUSTERS_FILE should look like is provided [here](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/clusters.tsv)


### Defragmentation

Most cloud genes in the pangenome are fragments of 'shell' or 'persistent' genes. Therefore, they do not provide informative data on the pangenome's diversity. 
To address this, we implemented an additional workflow to reduce the number of gene families and computational load by associating fragments with their original gene families.
This step is added to the previously described clustering process. 
The new workflow compares all representative protein sequences of gene families using the same identity threshold as the first step. 
It also uses the same coverage threshold, but only the smallest of the two protein sequences must be covered by at least the value specified by `--coverage`.

We then build a similarity graph, where the edges are the hits given by the comparison, and the nodes are the original gene families. 
Next, we iterate over all the nodes and compare them to their neighbors. 
If a node's neighbor has more members in its cluster and a longer representative sequence, we associate that node (and all the associated genes) with the neighbor. 
The genes linked to this node are considered as 'fragments' of the gene family represented by the neighboring node that is longer and more abundant.

To avoid using it, you can run the following:

```
ppanggolin cluster -p pangenome.h5 --no_defrag
```

In all cases, whichever pipeline you use, the gene families will end up in the 'pangenome.h5' file you enter as input.

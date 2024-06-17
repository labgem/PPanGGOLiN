### Cluster genes into gene families
 
After annotating genomes, their genes will be compared to determine similarities and build gene families using this information.

If .fasta files or annotation files containing gene sequences were provided, clustering can be executed by using the generated .h5 file, as such:

```
ppanggolin cluster -p pangenome.h5
```

#### How to customize MMSeqs2 clustering
```{warning}
Not all MMSeqs2 options are available in PPanGGOLiN. 
For a comprehensive overview of MMSeqs2 options, please refer to their documentation. 
To provide your own custom clustering from MMSeqs2 or another tool, please follow the instructions detailed in the [dedicated section](#read-clustering).
```

PPanGGOLiN will run [MMseqs2](https://github.com/soedinglab/MMseqs2) to perform clustering on all the protein sequences by searching for connected components for the clustering step.

##### How to set the identity and coverage parameters

PPanGGOLiN enables the setting of two essential parameters for gene clustering: **identity** and **coverage**. These parameters can be easily adjusted using `--identity` (default: 0.8) and `--coverage` (default: 0.8). The default values were selected as they are empirically effective parameters for aligning and clustering sequences at the species level.
 
Be aware that if you decrease identity and/or coverage, more genes will be clustered together in the same family. 
This will ultimately decrease the number of families and affect all subsequent steps.

```{note}
The chosen coverage mode in PPanGGOLiN requires both protein sequences to be covered by at least the proportion specified by --coverage, though this is modified afterwards by the [defragmentation step](./pangenomeAnalyses.md#defragmentation).
```

##### How to set the clustering mode

MMSeqs provides 3 different [clustering modes](https://github.com/soedinglab/MMseqs2/wiki#clustering-modes).
By default, the clustering mode used is the _single linkage_ (or _connected component_) algorithm.

Another option is the _set cover_ algorithm, which can be employed using `--mode 1`.

Additionally, the clustering algorithms of MMseqs, similar to CD-Hit, 
can be selected with `--mode 2` or its low-memory version through `--mode 3`.

(read-clustering)=
### Providing your gene families
 
It's possible to provide your clusters (or gene families); you must have provided the annotations in the first step. 
For gff3 files, the expected gene id is the 'ID' field in the 9<sup>th</sup> column. 
In the case of gbff or gbk files, use 'locus_tag' as a gene id, unless you are working with files from MaGe/MicroScope or SEED, where the id in the 'db_xref' field is used instead.

You must provide a .tsv file with a single gene id per line. This file must be at least with 2 columns.
The first column should indicate the cluster identifier or, as we call it, the gene family name.
The second column must correspond to the unique gene identifier used in the annotation files.
By default PPanGGOLiN will consider that the first line of a cluster (first occurrence of a family name) is also the line with the representative gene indicated. 
You can add a column with the representative gene identifier to use another representative gene.
You can also add a column to indicate if the gene is fragmented, with a 'F', or not (leave field empty).
*The column order is important !* You must first provide the cluster identifier, the representative id, and then the gene id to finish with the fragmented status of the gene.  

You can do this from the command line: 

`ppanggolin cluster -p pangenome.h5 --clusters clusters.tsv`

An example of what clusters.tsv should look like is provided [here](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/clusters.tsv)

```{note}
When you provide your clustering, *PPanGGOLiN* will translate the representative gene sequence of each family and write it in the HDF5 file.
```



### Defragmentation

Without performing additional steps, most cloud genes in the pangenome are fragments of 'shell' or 'persistent' genes. Therefore, they do not provide informative data on the pangenome's diversity. 
To address this, we implemented an additional step to the clustering to reduce the number of gene families and computational load by associating fragments to their original gene families.
This step is added to the previously described clustering process by default. 
It compares all representative protein sequences of gene families using the same identity threshold as the one given to MMseqs2, through  `--identity`. 
It also uses the same coverage threshold, but only the smallest of the two protein sequences must be covered by at least the value specified by `--coverage`.

We then build a similarity graph, where the edges are the hits given by the comparison, and the nodes are the original gene families. 
Next, we iterate over all the nodes and compare them to their neighbors. 
If a node's neighbor has more members in its cluster and a longer representative sequence, we associate that node (and all the associated genes) with the longer, more numerous neighbor. 
The genes linked to this node are considered as 'fragments' of the longer and more populated gene family represented by the neighboring node.

To avoid using this step, you can run the clustering with the following:
```
ppanggolin cluster -p pangenome.h5 --no_defrag
```
In all cases, whichever pipeline you use, the gene families will end up in the 'pangenome.h5' file you entered as input.
```{note}
this is performed only when running the clustering with PPanGGOLiN, and not when providing your own clustering results.
```

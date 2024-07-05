### Cluster genes into gene families

After annotating genomes, their genes will be compared to determine similarities and build gene families using this
information.

If .fasta files or annotation files containing gene sequences were provided, clustering can be executed by using the
generated .h5 file, as such:

```
ppanggolin cluster -p pangenome.h5
```

#### How to customize MMSeqs2 clustering

```{warning}
Not all MMSeqs2 options are available in PPanGGOLiN. 
For a comprehensive overview of MMSeqs2 options, please refer to their documentation. 
To provide your own custom clustering from MMSeqs2 or another tool, please follow the instructions detailed in the [dedicated section](#read-clustering).
```

PPanGGOLiN will run [MMseqs2](https://github.com/soedinglab/MMseqs2) to perform clustering on all the protein sequences
by searching for connected components for the clustering step.

##### How to set the identity and coverage parameters

PPanGGOLiN enables the setting of two essential parameters for gene clustering: **identity** and **coverage**. These
parameters can be easily adjusted using `--identity` (default: 0.8) and `--coverage` (default: 0.8). The default values
were selected as they are empirically effective parameters for aligning and clustering sequences at the species level.

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

It's possible to provide your own clusters (or gene families); you must have provided the annotations in the first step.
For gff3 files, the expected gene identifier is the 'ID' field in the 9<sup>th</sup> column. In the case of gbff or gbk
files, use 'locus_tag' as a gene identifier unless you are working with files from MaGe/MicroScope or SEED, where the id
in the 'db_xref' field is used instead.

You can give your clustering result in TSV file with a single gene identifier per line. This file consist of 2 to 4
columns in function of the information you want to give. See next part to have more information about expected format.


You can do this from the command line:

`ppanggolin cluster -p pangenome.h5 --clusters clusters.tsv`


If one of the gene in the pangenome is missing in your clustering, PPanGGOLiN will raise an error. To assert the gene into its own cluster (singleton) you can use the `--infer_singleton` option as such:

`ppanggolin cluster -p pangenome.h5 --clusters clusters.tsv --infer_singleton`

```{note}
When you provide your clustering, *PPanGGOLiN* will translate the representative gene sequence of each family and write it in the HDF5 file.
```


#### Assume the representative gene

The minimum information required is the gene family name (first column) and the constituent gene of the family (second
column). In that case, PPanGGOLiN will consider that the first line of a cluster (first occurrence of a family name) is
also the line with the representative gene indicated.

Here is a minimal example of your clustering file:

```
Family_A    Gene_1
Family_A    Gene_2
Family_A    Gene_3
Family_B    Gene_4
Family_B    Gene_5
Family_C    Gene_6
```

```{mermaid}

---
title: "Pangenome gene families when assuming representative gene"
align: center
---

%%{init: {'theme':'default'}}%%


flowchart TD
    subgraph C
        direction TB
        C6[6]
        style C6 fill:#00758f
    end
    style C rx:150,ry:150
    subgraph B
        direction TB
        B4[4]
        B5[5]
        style B4 fill:#00758f
    end
    style B rx:150,ry:150
    subgraph A
        direction TB
        A1[1]
        A2[2]
        A3[3]
        style A1 fill:#00758f
    end
    style A rx:150,ry:150
```

#### Specify the representative gene

It's possible to indicate which gene is the representative gene by adding a column between the gene family name.


Here is a minimal example of your clustering file:

```
Family_A    Gene_2  Gene_1
Family_A    Gene_2  Gene_2
Family_A    Gene_2  Gene_3
Family_B    Gene_5  Gene_4
Family_B    Gene_5  Gene_5
Family_C    Gene_6  Gene_6
```

```{mermaid}

---
title: "Pangenome gene families when assuming representative gene"
align: center
---

%%{init: {'theme':'default'}}%%


flowchart TD
    subgraph C
        direction TB
        C6[6]
        style C6 fill:#00758f
    end
    style C rx:150,ry:150
    subgraph B
        direction TB
        B4[4]
        B5[5]
        style B5 fill:#00758f
    end
    style B rx:150,ry:150
    subgraph A
        direction TB
        A1[1]
        A2[2]
        A3[3]
        style A2 fill:#00758f
    end
    style A rx:150,ry:150
```

#### Indicate fragmented gene

It's also possible to indicate if the gene is fragmented, by adding a new column in last position. Fragmented gene are tag with an 'F' in the last column.

You can add this column when you assume or not the representative gene. PPanGGOLiN will guess that this column is to precise the fragmented gene and assume if it must assert the representative gene

Here is a minimal example of your clustering file with fragmented gene precise:

```
Family_A    Gene_2  Gene_1  
Family_A    Gene_2  Gene_2  
Family_A    Gene_2  Gene_3  'F'
Family_B    Gene_5  Gene_4  
Family_B    Gene_5  Gene_5  
Family_C    Gene_6  Gene_6  'F'
```

```{warning
*The column order is
important !*  You must first provide the cluster identifier, the representative id, and then the gene id to finish with
the fragmented status of the gene.
```

### Defragmentation

Without performing additional steps, most cloud genes in the pangenome are fragments of 'shell' or 'persistent' genes.
Therefore, they do not provide informative data on the pangenome's diversity.
To address this, we implemented an additional step to the clustering to reduce the number of gene families and
computational load by associating fragments to their original gene families.
This step is added to the previously described clustering process by default.
It compares all representative protein sequences of gene families using the same identity threshold as the one given to
MMseqs2, through  `--identity`.
It also uses the same coverage threshold, but only the smallest of the two protein sequences must be covered by at least
the value specified by `--coverage`.

We then build a similarity graph, where the edges are the hits given by the comparison, and the nodes are the original
gene families.
Next, we iterate over all the nodes and compare them to their neighbors.
If a node's neighbor has more members in its cluster and a longer representative sequence, we associate that node (and
all the associated genes) with the longer, more numerous neighbor.
The genes linked to this node are considered as 'fragments' of the longer and more populated gene family represented by
the neighboring node.

To avoid using this step, you can run the clustering with the following:

```
ppanggolin cluster -p pangenome.h5 --no_defrag
```

In all cases, whichever pipeline you use, the gene families will end up in the 'pangenome.h5' file you entered as input.

```{note}
this is performed only when running the clustering with PPanGGOLiN, and not when providing your own clustering results.
```

(step-by-step-section)=
# Pangenome graph building and partitionning
The workflow subcommand of PPanGGOLiN uses what you can use through the other subcommands. In this section, the different steps will be briefly described along with the inputs that can be given, and which subcommands to use to run those specific parts yourself.

This is useful only if you want to customize your workflow parameters. If you want to build the pangenome of your species without tuning parameters, you can use the subcommand `ppanggolin workflow` as described in the introduction.

## Genomes annotation and storage

### Annotate pangenome with fasta files

As an input file, you can provide a list of .fasta files. If you do so, the provided genomes will be annotated using the following tools:  The CDS will be annotated using [Prodigal](https://github.com/hyattpd/Prodigal), the tRNA will be annotated using [ARAGORN](http://130.235.244.92/ARAGORN/) and the rRNA are annotated using the [Infernal](http://eddylab.org/infernal/) command-line tools coupled with HMM of the bacterial and archaeal rRNAs downloaded from [RFAM](https://rfam.xfam.org/). Then the CDS overlapping any RNA genes will be deleted as they are usually false positive calls.

To run this part of the pipeline, you can do : 

`ppanggolin annotate --fasta ORGANISM_FASTA_LIST`

With ORGANISM_FASTA_LIST following the format described for the 'workflow' subcommand. You can check [this example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.fasta.list).

### Annotate pangenome with annotation files

If you do not want to use this pipeline, you can provide your annotation files (This is especially recommended if you already have functional annotations of your genome, as they will be added to the pangenome).

They can be either gff3 files or .gbk files or .gbff files, or a mix of them, and should be provided through a list alike [this example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list). .gbk or .gbff files are preferred.

You can provide them using the following command : 

`ppanggolin annotate --anno ORGANISM_ANNOTATION_LIST`

With ORGANISM_ANNOTATION_LIST being your file listing the organisms and the annotation file associated. If your annotation files do not have the genome sequence in them, you can use both options at the same time (to have both the gene annotations and the gene sequences) as such : 

`ppanggolin annotate --anno ORGANISM_ANNOTATION_LIST --fasta ORGANISM_FASTA_LIST`

In addition, you can tune the command that is run with the annotation pipeline, or read in your annotation files using various options. You can check them through the command line `ppanggolin annotate --help`, they help should be self-explanatory. If not, don't hesitate to ask questions through the [issues page](https://github.com/labgem/PPanGGOLiN/issues).

### Annotate command-line options

It's possible to tune the annotation command with some parameters describe as follows: 
```bash
  --allow_overlap       Use to not remove genes overlapping with RNA features. (default: False)
  --norna               Use to avoid annotating RNA features. (default: False)
  --kingdom {bacteria,archaea}
                        Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. (default: bacteria)
  --translation_table TRANSLATION_TABLE
                        Translation table (genetic code) to use. (default: 11)
  --basename BASENAME   basename for the output file (default: pangenome)
  --use_pseudo          In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them)
                        (default: False)
  -p {single,meta}, --prodigal_procedure {single,meta}
                        Allow to force the prodigal procedure. If nothing given, PPanGGOLiN will decide in function of contig length (default:
                        None)
```


## Clustering genes in gene families

### Computing gene families with PPanGGOLiN

Once we have the genes, we need to compare them to know which are similar, and to build gene families through this information. 
If you provided .fasta files or annotation files with gene sequences in them, clustering can be run directly by providing the .h5 file that was generated, as such : 

`ppanggolin cluster -p pangenome.h5`

PPanGGOLiN will call [MMseqs2](https://github.com/soedinglab/MMseqs2) to run the clustering on all the protein sequences by searching for connected components for the clustering step. You can tune its parameters using `--identity`(default 0.8) and `--coverage`(default 0.8). You can use other clustering algorithms of MMseqs by using --mode (default 1). Both protein sequences have to be covered by at least the proportion indicated by --coverage.

### Providing your gene families

If you do not want to use MMseqs2 and provide your clusters (or gene families) you can do so only if you provided the annotations in the first step. In the case of gff3 files, the 'ID' field in the 9th column is expected as a gene id. In the case of gbff or gbk files, the 'locus_tag' is used as a gene id, except with files coming from MaGe or from SEED, where the id provided in the 'db_xref' field is used.

You will need to provide a .tsv file. The first column indicates the cluster id, and the second column indicates a unique gene id that is used in the annotation files. There is a single gene id per line.

You can do that through the command line : 

`ppanggolin cluster -p pangenome.h5 --clusters MY_CLUSTERS_FILE`

An example of what MY_CLUSTERS_FILE should look like is provided [here](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/clusters.tsv)

There are other options that you can use to tune your clustering. Most of them should be self-explanatory. If not, do not hesitate to write an [issue](https://github.com/labgem/PPanGGOLiN/issues). The only tricky option is the '--no-defrag' option. 

### Defragmentation

We noticed that most of the cloud genes in the pangenome are fragments of 'shell' or 'persistent' genes, and so not informative on the pangenome's diversity. We added an additional workflow to reduce the number of gene families and reduce the computational load by trying to associate fragments to their original gene families.

It adds a step to the clustering described previously. It will compare all of the gene families representative protein sequences using the same identity threshold as the first step. It will also use the same coverage threshold, but only the smallest of both protein sequences have to be covered by at least the value indicated by '--coverage'.

After that, we build a similarity graph where the edges are the hits given by the comparison, and the nodes are the original gene families. Then we iter on all nodes and compare them to their neighbors. If the neighbor of a node is more numerous (has more members in the cluster it represents) and its representative sequence is longer, that node (and all the genes associated) is associated with the neighbor. The genes associated with this node are defined as 'fragments' of the gene family represented by the longer and more numerous neighboring node.

You can use this pipeline by using the following command for versions 1.1.88 and under: 

`ppanggolin cluster -p pangenome.h5 --defrag`

Starting from 1.1.89 and on, this strategy has become the default behavior of PPanGGOLiN. To avoid using it, you can run the following:

`ppanggolin cluster -p pangenome.h5 --no_defrag`

In any case and whichever pipeline you use, in the end, the gene families will be saved in the 'pangenome.h5' given as input.

## Building the pangenome graph

To partition a pangenome graph, you need to build a said pangenome graph. This can be done through the 'graph' subcommand. This will take a pangenome .h5 file as input and add edges to it.
This subcommand has only a single other option, which is '-r' or '--remove_high_copy_number'. If used, it will remove the gene families that are too duplicated in your genomes. This is useful if you want to visualize your pangenome afterward and want to remove the biggest hubs to have a clearer view. It can also be used to limit the influence of very duplicated genes such as transposase or ABC transporters in the partition step.

The graph is constructed using the following subcommand : 

`ppanggolin graph -p pangenome.h5`

And the results are saved in the pangenome.h5 file given as input.

## Pangenome graph partitioning 

This is the step that will assign gene families to the 'persistent', 'shell', or 'cloud' partitions. 

The 'persistent' partition will group genes that are present throughout the entire species. They will be essential genes, genes required for important metabolic pathways and genes that define the metabolic and biosynthetic capabilities of the taxonomic group.

The 'shell' partition groups genes that are present in only some individuals. Those are often genes that were acquired through horizontal gene transfers and encode for functions involved in environmental adaptations, pathogenicity, virulence or encoding secondary metabolites for example.

The 'cloud' partition groups genes that are very rare in the pangenome and found in one, or very few, individuals. Most of the genes were associated with phage-related genes. They probably all were acquired through horizontal gene transfers. Antibiotic resistance genes were often found to be belonging to the cloud genome, as well as plasmid genes.

It can be realized through the following subcommand : 

`ppanggolin partition -p pangenome.h5`

It also has quite a few options. Most of them are not self-explanatory. If you want to know what they do, you should read the PPanGGOLiN paper (you can read it [here](https://journals.plos.org/ploscompbiol/article?rev=2&id=10.1371/journal.pcbi.1007732)) where the statistical methods used are thoroughly described.

The one parameter that might be of importance is the '-K', or '--nb_of_partitions' parameter. This will define the number of classes used to partition the pangenome. This may be of use if you expect to have well-defined subpopulations in your pangenome, and you know exactly how many. If not, that number is detected automatically through an ICL criterion. The idea is that the most present partition will be 'persistent', the least present will be 'cloud', and all the others will be 'shell'. The number of partitions corresponding to the shell will be the number of expected subpopulations in your pangenome. (So if you expect 5 subpopulations, you could use -K 7). 

In most cases, you should let the statistical criterion used by PPanGGOLiN find the optimal number of partitions for you.

All the results will be added to the given 'pangenome.h5' input file.
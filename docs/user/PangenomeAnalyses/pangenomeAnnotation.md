(annot-fasta)=
### Annotate fasta file

As an input file, you can provide a list of .fasta files. 
If you do so, the provided genomes will be annotated using the following tools: 

- [Pyrodigal](https://pyrodigal.readthedocs.io/en/stable/index.html) , which is based on Prodigal, to annotate CDSs
- [ARAGORN](http://www.ansikte.se/ARAGORN/) to annotate  tRNAs
- [Infernal](http://eddylab.org/infernal/) coupled with HMM of the bacterial and archaeal rRNAs downloaded from [RFAM](https://rfam.xfam.org/) to annotate RNAs.

To proceed with this stage of the pipeline, you need to create an ORGANISMS_FASTA_LIST. This file should be tab-separated with each line depicting an individual organism and its pertinent information with the following organisation:

- The first column contains a unique organism name
- The second column contains the path to the associated FASTA file and circular
- The following columns contain Contig identifiers 
   

You can check [this example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.fasta.list)

To run the annotation part, you can use this minimal command:

```
ppanggolin annotate --fasta ORGANISM_FASTA_LIST
```

#### Use a different genetic code in my annotation step
To annotate the genomes, you can easily change the translation table (or genetic code) used by Pyrodigal just by giving the corresponding number as described [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

#### Force the Prodigal procedure
Prodigal can predict gene in [single/normal mode](https://github.com/hyattpd/prodigal/wiki/gene-prediction-modes#normal-mode) that include a training step on your genomes or in [meta/anonymous mode](https://github.com/hyattpd/prodigal/wiki/gene-prediction-modes#anonymous-mode) which use pre-calculated training files. 
As recommended in the Prodigal documentation: "Anonymous mode should be used on metagenomic data sets, or on sequences too short to provide good training data."

By default PPanGGOLiN will determine the best mode based on the contig length.
The procedure can be overriden with the option `-p, --prodigal_procedure`.
The option only accepts  **single** or **meta** keywords, corresponding to the Prodigal procedure name.

#### Customize the RNA annotation
If you do not want to predict the RNA (and thus not use Infernal and Aragorn), you can add the `--norna` option to your command.
Otherwise, by default, any CDS overlapping  RNA genes will be deleted as they are often false positive calls.
You can prevent this filtering by using the `--allow_overlap` option.

Additionally, the `--kingdom archaea` option can be utilized when working with archaea genomes to specify infernal's RNA annotation model. 

### Use annotated file as pangenome base

You can  provide  annotation files in various formats inclusding gff3 files or .gbk files or .gbff files, or a mix of them. They should be provided through as a list in a tab-separated file similar to [this example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list).

```{note}
Use your own annotation for your genome is highly recommended, particularly if you already
 have functional annotations, as they can be added to the pangenome.
```

You can provide them using the following command: 

```
ppanggolin annotate --anno ORGANISM_ANNOTATION_LIST
```

#### How to deal with annotation files without sequences

If your annotation files do  not contain the genome sequence, 
you can use the both options simultaneously to obtain the gene annotations and gene sequences, as follows: 

```
ppanggolin annotate --anno ORGANISM_ANNOTATION_LIST --fasta ORGANISM_FASTA_LIST
```

#### Take the pseudogenes into account for pangenome analyses

By default PPanGGOLiN will not take pseudogenes into account. However, they could be interesting in some context.
It is possible to include pseudogenes in the pangenome by using the `--use_pseudo`option.

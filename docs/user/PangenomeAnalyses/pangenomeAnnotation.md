(annot-fasta)=
### Annotate fasta file

As an input file, you can provide a list of .fasta files. 
If you do so, the provided genomes will be annotated using the following tools: 

- [Pyrodigal](https://pyrodigal.readthedocs.io/en/stable/index.html) to annotate the CDS, which is based on Prodigal,
- [ARAGORN](http://www.ansikte.se/ARAGORN/) to annotate the tRNA
- [Infernal](http://eddylab.org/infernal/) coupled with HMM of the bacterial and archaeal rRNAs downloaded from [RFAM](https://rfam.xfam.org/) to annotate the RNA command-line tools.

To run this part of the pipeline, you must create an ORGANISMS_FASTA_LIST which is a tab-separated file with the following organisation :

1. The first column contains a unique organism name
2. The second column the path to the associated FASTA file
3. Circular contig identifiers are indicated in the following columns
4. Each line represents an organism

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
By default PPanGGOLiN will decide the best mode in function of the contig length.
It's possible to force the procedure with the option `-p, --prodigal_procedure`.
The option accept only **single** or **meta** keyword, corresponding to the prodigal procedure name.

#### Customize the RNA annotation
If you don't want to predict the RNA (and so don't use Infernal and Aragorn) you can add the option `--norna` in your command.
Else, by default the CDS overlapping any RNA genes will be deleted as they are usually false positive calls.
You can prevent this filtering by using `--allow_overlap` option.

Moreover, when you are working with archaea genomes, you can use the option `--kingdom archaea` to indicate to infernal which model to use to annotate RNA. 

### Use annotated file as pangenome base

You can also provide your annotation files.
They can be either gff3 files or .gbk files or .gbff files, or a mix of them, and should be provided through a list in a tab-separated file alike [this example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list).

```{note}
Use your own annotation is especially recommended if you already have functional annotations of your genome,
 as they will be added to the pangenome
```

You can provide them using the following command : 

```
ppanggolin annotate --anno ORGANISM_ANNOTATION_LIST
```

#### How to deal with annotation files without sequences

If your annotation files do not have the genome sequence in them, 
you can use both options at the same time (to have both the gene annotations and the gene sequences) as such : 

```
ppanggolin annotate --anno ORGANISM_ANNOTATION_LIST --fasta ORGANISM_FASTA_LIST
```

#### Take the pseudogenes into account for pangenome analyses

By default PPanGGOLiN will not take into account the pseudogene. However, they could be interesting in some context.
So it's possible to add the pseudogenes in the pangenome with the option `--use_pseudo`.

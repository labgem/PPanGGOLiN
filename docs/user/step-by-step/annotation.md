### Annotate pangenome with fasta files

As an input file, you can provide a list of .fasta files. 
If you do so, the provided genomes will be annotated using the following tools: 

- [Prodigal](https://github.com/hyattpd/Prodigal) to annotate the CDS,
- [ARAGORN](http://130.235.244.92/ARAGORN/) to annotate the tRNA
- [Infernal](http://eddylab.org/infernal/) coupled with HMM of the bacterial and archaeal rRNAs downloaded from [RFAM](https://rfam.xfam.org/) to annotate the RNA command-line tools.

Then the CDS overlapping any RNA genes will be deleted as they are usually false positive calls.
You can prevent this filtering by using `--allow_overlap` option.

To run this part of the pipeline, you can do :

```
ppanggolin annotate --fasta ORGANISM_FASTA_LIST
```

The file ORGANISMS_FASTA_LIST is a tsv-separated file with the following organisation :

1. The first column contains a unique organism name
2. The second column the path to the associated FASTA file
3. Circular contig identifiers are indicated in the following columns
4. Each line represents an organism

You can check [this example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.fasta.list).

### Annotate pangenome with annotation files

You can also provide your annotation files.
They can be either gff3 files or .gbk files or .gbff files, or a mix of them, and should be provided through a list alike [this example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list).
.gbk or .gbff files are preferred.

```{note}
Use your own annotation is especially recommended if you already have functional annotations of your genome,
 as they will be added to the pangenome
```

You can provide them using the following command : 

```
ppanggolin annotate --anno ORGANISM_ANNOTATION_LIST
```

With ORGANISM_ANNOTATION_LIST being your file listing the organisms and the annotation file associated. 
If your annotation files do not have the genome sequence in them, 
you can use both options at the same time (to have both the gene annotations and the gene sequences) as such : 

```
ppanggolin annotate --anno ORGANISM_ANNOTATION_LIST --fasta ORGANISM_FASTA_LIST
```

### Annotate command-line options

You can tune the command that is run with the annotation pipeline, or read in your annotation files using various options described below.

| name                 | alias | default   | type / choices     | description                                                                                                     |
|----------------------|-------|-----------|--------------------|-----------------------------------------------------------------------------------------------------------------|
| --allow_overlap      |       | False     | bool               | Use to not remove genes overlapping with RNA features                                                           | 
| --norna              |       | False     | bool               | Use to avoid annotating RNA features                                                                            |
| --kingdom            |       | bacteria  | {bacteria,archaea} | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation                     |
| --translation_table  |       | 11        | integer            | Translation table (genetic code) to use                                                                         |
| --basename           |       | pangenome | string             | basename for the output file                                                                                    |
| --use_pseudo         |       | False     | bool               | In the context of provided annotation, use this option to read pseudogenes (Default behavior is to ignore them) |
| --prodigal_procedure | -p    | None      | {single,meta}      | Allow to force the prodigal procedure. If nothing given, PPanGGOLiN will decide in function of contig length    |
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
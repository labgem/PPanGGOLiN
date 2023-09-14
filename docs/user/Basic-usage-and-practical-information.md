(basic)=
# Basic usage and pratical information

## The 'workflow' subcommand

We tried to make PPanGGOLiN relatively easy to use by making this **'workflow'** subcommand. 
It runs a pangenome analysis whose exact steps will depend on the input files you provide it with. 
In the end, you will end up with some files and figures that describe the pangenome of your taxonomic group of interest in different ways.

The minimal subcommand is as follows :
 
```
ppanggolin workflow --fasta ORGANISMS_FASTA_LIST
```

It uses parameters that we found to be generally the best when working with species pangenomes.

The file ORGANISMS_FASTA_LIST is a tsv-separated file with the following organisation :

1. The first column contains a unique organism name
2. The second column the path to the associated FASTA file
3. Circular contig identifiers are indicated in the following columns
4. Each line represents an organism


An [example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.fasta.list) with 50 *Chlamydia trachomatis* genomes can be found in the testingDataset/ directory.

You can also give PPanGGOLiN your own annotations using .gff or .gbff/.gbk files instead of .fasta files as long as they include the genomic dna sequences, such as the ones provided by prokka using the following command :

```
ppanggolin workflow --anno ORGANISMS_ANNOTATION_LIST
```

Another [example](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/organisms.gbff.list) of such a file can be found in the testingDataset/ directory.


```{note}
Look at the **annotate** command documentation for more information [here](#annotation)
```

In addition, you can provide your own gene families. 
PPanGGOLiN will use it to build and partition the pangenome graph.
You can do that through the command line : 

```
ppanggolin workflow --fasta ORGANISMS_FASTA_LIST --anno ORGANISMS_ANNOTATION_LIST --clusters MY_CLUSTERS_FILE
```

An example of what MY_CLUSTERS_FILE should look like is provided [here](https://github.com/labgem/PPanGGOLiN/blob/master/testingDataset/clusters.tsv)

Whether you use fasta or annotations, the workflow command options are the same.

| name              | alias | default                         | type / choices     | description                                                                                                                   |
|-------------------|-------|---------------------------------|--------------------|-------------------------------------------------------------------------------------------------------------------------------|
| output            | -o    | ppanggolin_output_DATE_HOUR_PID | Path               | Output directory to save the pangenome and all the output files                                                               |
| basename          |       | pangenome                       | string             | basename for the pangenome file                                                                                               |
| rarefaction       |       | False                           | bool               | Use to compute the rarefaction curves (WARNING: can be time consuming)                                                        |
| cpu               | -c    | 1                               | integer            | Number of available cpus                                                                                                      |
| translation_table |       | 11                              | integer            | Translation table (genetic code) to use                                                                                       |
| kingdom           |       | bacteria                        | {bacteria,archaea} | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation                                   |
| mode              |       | 1                               | {0,1,2,3}          | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) |
| coverage          |       | 0.8                             | 0<=float<=1        | Minimal coverage of the alignment for two proteins to be in the same cluster                                                  |
| identity          |       | 0.8                             | 0<=float<=1        | Minimal identity percent for two proteins to be in the same cluster                                                           |
| nb_of_partitions  | -K    | -1                              | integer            | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically                                |
| no_defrag         |       | False                           | bool               | DO NOT Realign gene families to link fragments with their non-fragmented gene family                                          |
| no_flat_files     |       | False                           | bool               | Generate only the HDF5 pangenome file                                                                                         |
| tmpdir            |       | TMPDIR                          | Path               | directory for storing temporary files                                                                                         |

(panrgp)=
## The 'panrgp' subcommand 

This command works exactly like 'workflow'. The difference is that it will run more analysis related to [Regions of Genome Plasticity](#RGP-section).
You can use the panrgp command as follow:

```bash
ppanggolin panrgp --fasta ORGANISMS_FASTA_LIST
```

The rgp analysis is launched after the pangenome partitionning and use the default parameters. 
If you want to tune the rgp detection, you can use the `rgp` command after the `workflow` command.


More detail about RGP detection [here](#RGP-section) and in the [panRGP publication](https://doi.org/10.1093/bioinformatics/btaa792)

(panmodule)=
## The 'panmodule' subcommand

Again, it works like 'workflow' but you can detect the conserved modules in your pangenome, you can use the **panModule** workflow, as such:

```bash
ppanggolin panmodule --fasta ORGANISMS_FASTA_LIST
```

The module prediction is launched after the pangenome partitionning with the default parameters. 
If you want to tune the module detection, you can use the `module` command after the `workflow`.


Further details can be found in the [conserved module analysis documentation](#module-section) and in the [panModule publication](https://doi.org/10.1101/2021.12.06.471380)

## Run all PPanGGOLiN analysis

Finally it's also possible to run all analysis with one command wrapper `all`. 
With this workflow, the pangenome will be built and partionned and RGP, spots and module will be predicted. 
You can run all the analysis as such:

```bash
ppanggolin all --fasta ORGANISMS_FASTA_LIST
```

## Configuration file

Advanced users can provide a configuration file containing any or all parameters to PPanGGolin commands. 
This feature is particularly useful for workflow commands such as `workflow`, `all`, `panrgp`, and `panmodule`, as it allows for the specification of all parameters for each subcommand launched in a workflow. 
Additionally, a configuration file can be used to reuse a specific set of parameters across multiple pangenomes.

To provide a configuration file to a PPanGGolin command, use the `--config` parameter.

```{note} 
Any command line arguments provided along with a configuration file will override the corresponding arguments specified in the configuration file.
When an argument is not specified in either the command line or the configuration file, the default value is used.
```

The configuration file is a JSON file that contains two sections common to all commands: `input_parameters` and `general_parameters`. 
In addition, there is a section for each subcommand that contains its specific parameters.

You can generate a configuration file template with default values by using the `ppanggolin utils` command as follows:

```
ppanggolin utils --default_config CMD
```

For example, to generate a configuration file for the panrgp command with default values, use the command 
```
ppanggolin utils --default_config panrgp
```
 
 This command will create the following configuration file: 

```yaml
input_parameters:
    # A tab-separated file listing the organism names, and the fasta filepath of its
    # genomic sequence(s) (the fastas can be compressed with gzip). One line per organism.
  # fasta: <fasta file>
    # A tab-separated file listing the organism names, and the gff/gbff filepath of
    # its annotations (the files can be compressed with gzip). One line
    # per organism. If this is provided, those annotations will be used.
  # anno: <anno file>

general_parameters:
    # Output directory
    output: ppanggolin_output_DATE2023-04-14_HOUR10.09.27_PID14968
    # basename for the output file
    basename: pangenome
    # directory for storing temporary files
    tmpdir: /tmp
    # Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)
    # Choices: 0, 1, 2
    verbose: 1
    # log output file
    log: stdout
    # disables the progress bars
    disable_prog_bar: False
    # Force writing in output directory and in pangenome output file.
    force: False

annotate:
    # Use to not remove genes overlapping with RNA features.
    allow_overlap: False
    # Use to avoid annotating RNA features.
    norna: False
    # Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.
    # Choices: bacteria, archaea
    kingdom: bacteria
    # Translation table (genetic code) to use.
    translation_table: 11
    # In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them)
    use_pseudo: False
    # Allow to force the prodigal procedure. If nothing given, PPanGGOLiN will decide in function of contig length
    # Choices: single, meta
    prodigal_procedure: False
    # Number of available cpus
    cpu: 1
```

## Required computing resources

Most of PPanGGOLiN's commands should be run with as many CPUs as you can give them by using the --cpu option as PPanGGOLiN's speed increases relatively well with the number of CPUs. 
While the 'smallest' pangenomes (up to a few hundred genomes) can be easily analyzed on a normal desktop computer, 
the biggest ones will require a good amount of RAM.
For example, 40 strains of *E. coli* were analyzed in 3 minutes using 1.2Go of RAM using 16 threads. 
1000 strains were analyzed in 45 minutes with 14 Go of RAM using 16 threads, and as of writing those lines, 
20 656 genomes was the biggest pangenome we did, and it required about a day and 120 Go of RAM.
The following graphic can give you an idea of the time it takes for a pangenome analysis given the number of genomes in input:

```{image} ../_static/runtimes.png
:align: center
```

## Usage and basic options

As most programs in bioinformatics, you can always specify some utility options.

You can specify the number of CPUs to use (which is recommended ! The default is to use just one) using the option `--cpu`.

You can specify the output directory (if not provided, one will be generated) using the option `--output`.

If you work in a strange environment that has no, or little available disk space in the '/tmp' (or your system equivalent, what is stored in TMPDIR) directory, you can specify a new temporary directory using `--tmp`

If you want to redo an analysis from scratch and store it in a directory that already exists, you will have to use the `--force` option. 
Be wary, however, that the data in that directory will be overwritten if named identically as any output file written by ppanggolin.

PPanGGOLiN is deliberately very verbose, to help users understand each stage of the analysis. 
If you want, verbosity can be reduced in several ways.
First, you can specify the verbosity level with the `--verbose` option. 
With `0` will show only warning and erros, `1` will add the information (default value), and if you encounter any problem you can use the debug level with value `2`.
Then you can also remove the progress bar with the option `--disable_prog_bar`
Finaly, you can also save PPanGGOLiN logs in a file by specified its path with the option `--log`.
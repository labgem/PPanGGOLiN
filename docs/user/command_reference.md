# Command-line reference

Complete reference of the PPanGGOLiN command-line interface.

This page is automatically generated from the PPanGGOLiN command-line parser.

## Basic

### `ppanggolin all`

Easy workflow to run all possible analysis

#### Input arguments for ppanggolin all

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fasta` | Path | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per genome. This option can be used alone. |
| `--anno` | Path | — | A tab-separated file listing the genome names, and the gff filepath of its annotations (the gffs can be compressed). One line per genome. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used. |
| `--clusters` | Path | — | a tab-separated file listing the cluster names, the gene IDs, and optionally whether they are a fragment or not. |

#### Optional arguments for ppanggolin all

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `--basename` | str | `pangenome` | basename for the output file |
| `--rarefaction` | bool | False | Use to compute the rarefaction curves (WARNING: can be time consuming) |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--kingdom` | lower | `bacteria` | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. <br>Choices: `bacteria`, `archaea` |
| `--mode` | str | `1` | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) <br>Choices: `0`, `1`, `2`, `3` |
| `--coverage` | float | 0.8 | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--identity` | float | 0.8 | Minimal identity percent for two proteins to be in the same cluster |
| `--infer_singletons` | bool | False | Use this option together with --clusters. If a gene is not present in the provided clustering result file, it will be assigned to its own unique cluster as a singleton. |
| `--use_pseudo` | bool | False | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-K, --nb_of_partitions` | int | -1 | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `--no_defrag` | bool | False | DO NOT Realign gene families to link fragments with their non-fragmented gene family. |
| `--no_flat_files` | bool | False | Generate only the HDF5 pangenome file. |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |

#### Common arguments for ppanggolin all

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin panmodule`

Easy workflow to run a pangenome analysis with module prediction

#### Input arguments for ppanggolin panmodule

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fasta` | Path | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per genome. This option can be used alone. |
| `--anno` | Path | — | A tab-separated file listing the genome names, and the gff filepath of its annotations (the gffs can be compressed). One line per genome. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used. |
| `--clusters` | Path | — | a tab-separated file listing the cluster names, the gene IDs, and optionally whether they are a fragment or not. |

#### Optional arguments for ppanggolin panmodule

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `--basename` | str | `pangenome` | basename for the output file |
| `--rarefaction` | bool | False | Use to compute the rarefaction curves (WARNING: can be time consuming) |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--kingdom` | lower | `bacteria` | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. <br>Choices: `bacteria`, `archaea` |
| `--mode` | str | `1` | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) <br>Choices: `0`, `1`, `2`, `3` |
| `--coverage` | float | 0.8 | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--identity` | float | 0.8 | Minimal identity percent for two proteins to be in the same cluster |
| `--infer_singletons` | bool | False | Use this option together with --clusters. If a gene is not present in the provided clustering result file, it will be assigned to its own unique cluster as a singleton. |
| `--use_pseudo` | bool | False | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-K, --nb_of_partitions` | int | -1 | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `--no_defrag` | bool | False | DO NOT Realign gene families to link fragments with their non-fragmented gene family. |
| `--no_flat_files` | bool | False | Generate only the HDF5 pangenome file. |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |

#### Common arguments for ppanggolin panmodule

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin panrgp`

Easy workflow to run a pangenome analysis with genomic islands and spots of insertion detection

#### Input arguments for ppanggolin panrgp

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fasta` | Path | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per genome. This option can be used alone. |
| `--anno` | Path | — | A tab-separated file listing the genome names, and the gff filepath of its annotations (the gffs can be compressed). One line per genome. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used. |
| `--clusters` | Path | — | a tab-separated file listing the cluster names, the gene IDs, and optionally whether they are a fragment or not. |

#### Optional arguments for ppanggolin panrgp

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `--basename` | str | `pangenome` | basename for the output file |
| `--rarefaction` | bool | False | Use to compute the rarefaction curves (WARNING: can be time consuming) |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--kingdom` | lower | `bacteria` | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. <br>Choices: `bacteria`, `archaea` |
| `--mode` | str | `1` | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) <br>Choices: `0`, `1`, `2`, `3` |
| `--coverage` | float | 0.8 | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--identity` | float | 0.8 | Minimal identity percent for two proteins to be in the same cluster |
| `--infer_singletons` | bool | False | Use this option together with --clusters. If a gene is not present in the provided clustering result file, it will be assigned to its own unique cluster as a singleton. |
| `--use_pseudo` | bool | False | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-K, --nb_of_partitions` | int | -1 | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `--no_defrag` | bool | False | DO NOT Realign gene families to link fragments with their non-fragmented gene family. |
| `--no_flat_files` | bool | False | Generate only the HDF5 pangenome file. |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |

#### Common arguments for ppanggolin panrgp

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin workflow`

Easy workflow to run a pangenome analysis in one go

#### Input arguments for ppanggolin workflow

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fasta` | Path | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per genome. This option can be used alone. |
| `--anno` | Path | — | A tab-separated file listing the genome names, and the gff filepath of its annotations (the gffs can be compressed). One line per genome. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used. |
| `--clusters` | Path | — | a tab-separated file listing the cluster names, the gene IDs, and optionally whether they are a fragment or not. |

#### Optional arguments for ppanggolin workflow

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `--basename` | str | `pangenome` | basename for the output file |
| `--rarefaction` | bool | False | Use to compute the rarefaction curves (WARNING: can be time consuming) |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--kingdom` | lower | `bacteria` | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. <br>Choices: `bacteria`, `archaea` |
| `--mode` | str | `1` | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) <br>Choices: `0`, `1`, `2`, `3` |
| `--coverage` | float | 0.8 | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--identity` | float | 0.8 | Minimal identity percent for two proteins to be in the same cluster |
| `--infer_singletons` | bool | False | Use this option together with --clusters. If a gene is not present in the provided clustering result file, it will be assigned to its own unique cluster as a singleton. |
| `--use_pseudo` | bool | False | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-K, --nb_of_partitions` | int | -1 | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `--no_defrag` | bool | False | DO NOT Realign gene families to link fragments with their non-fragmented gene family. |
| `--no_flat_files` | bool | False | Generate only the HDF5 pangenome file. |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |

#### Common arguments for ppanggolin workflow

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


## Expert

### `ppanggolin annotate`

Annotate genomes

#### Required arguments for ppanggolin annotate

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fasta` | Path | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed with gzip). One line per genome. |
| `--anno` | Path | — | A tab-separated file listing the genome names, and the gff/gbff filepath of its annotations (the files can be compressed with gzip). One line per genome. If this is provided, those annotations will be used. |

#### Optional arguments for ppanggolin annotate

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `--allow_overlap` | bool | False | Use to not remove genes overlapping with RNA features. |
| `--norna` | bool | False | Use to avoid annotating RNA features. |
| `--kingdom` | lower | `bacteria` | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. <br>Choices: `bacteria`, `archaea` |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--basename` | str | `pangenome` | basename for the output file |
| `--use_pseudo` | bool | False | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-p, --prodigal_procedure` | lower | — | Allow to force the prodigal procedure. If nothing given, PPanGGOLiN will decide in function of contig length <br>Choices: `single`, `meta` |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |

#### Common arguments for ppanggolin annotate

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin cluster`

Cluster genes into gene families

#### Required arguments for ppanggolin cluster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |

#### Clustering arguments for ppanggolin cluster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--identity` | float | 0.8 | Minimal identity percent for two proteins to be in the same cluster |
| `--coverage` | float | 0.8 | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--mode` | str | `1` | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) <br>Choices: `0`, `1`, `2`, `3` |
| `--no_defrag` | bool | False | DO NOT Use the defragmentation strategy to link potential fragments with their original gene family. |

#### Read clustering arguments for ppanggolin cluster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--clusters` | Path | — | A tab-separated list containing the result of a clustering. One line per gene. First column is cluster ID, and second is gene ID |
| `--infer_singletons` | bool | False | When reading a clustering result with --clusters, if a gene is not in the provided file it will be placed in a cluster where the gene is the only member. |

#### Optional arguments for ppanggolin cluster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--tmpdir` | Path | `/tmp` | directory for storing temporary files |
| `--keep_tmp` | bool | False | Keeping temporary files (useful for debugging). |

#### Common arguments for ppanggolin cluster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin graph`

Create the pangenome graph

#### Required arguments for ppanggolin graph

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |

#### Optional arguments for ppanggolin graph

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-r, --remove_high_copy_number` | int | 0 | Positive Number: Remove families having a number of copy of gene in a single genome above or equal to this threshold in at least one genome (0 or negative values are ignored). |

#### Common arguments for ppanggolin graph

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin metadata`

Add metadata to elements in yout pangenome

#### Required arguments for ppanggolin metadata

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |
| `-m, --metadata` | Path | — | Metadata in TSV file. See our github for more detail about format |
| `-s, --source` | str | — | Name of the metadata source |
| `-a, --assign` | str | — | Select to which pangenome element metadata will be assigned <br>Choices: `families`, `genomes`, `contigs`, `genes`, `RGPs`, `spots`, `modules` |

#### Optional arguments for ppanggolin metadata

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--omit` | bool | False | Allow to pass if a key in metadata is not find in pangenome |

#### Common arguments for ppanggolin metadata

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin partition`

Partition the pangenome graph

#### Required arguments for ppanggolin partition

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome.h5 file |

#### Optional arguments for ppanggolin partition

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-b, --beta` | float | 2.5 | beta is the strength of the smoothing using the graph topology during partitioning. 0 will deactivate spatial smoothing. |
| `-ms, --max_degree_smoothing` | float | 10 | max. degree of the nodes to be included in the smoothing process. |
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `-fd, --free_dispersion` | bool | False | use if the dispersion around the centroid vector of each partition during must be free. It will be the same for all genomes by default. |
| `-ck, --chunk_size` | int | 500 | Size of the chunks when performing partitioning using chunks of genomes. Chunk partitioning will be used automatically if the number of genomes is above this number. |
| `-K, --nb_of_partitions` | int | -1 | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `-Kmm, --krange` | int | `[3, 20]` | Range of K values to test when detecting K automatically. |
| `-im, --ICL_margin` | float | 0.05 | K is detected automatically by maximizing ICL. However at some point the ICL reaches a plateau. Therefore we are looking for the minimal value of K without significant gain from the larger values of K measured by ICL. For that we take the lowest K that is found within a given 'margin' of the maximal ICL value. Basically, change this option only if you truly understand it, otherwise just leave it be. |
| `--draw_ICL` | bool | False | Use if you want to draw the ICL curve for all the tested K values. Will not be done if K is given. |
| `--keep_tmp_files` | bool | False | Use if you want to keep the temporary NEM files |
| `-se, --seed` | int | 42 | seed used to generate random numbers |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |

#### Common arguments for ppanggolin partition

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin rarefaction`

Compute the rarefaction curve of the pangenome

#### Required arguments for ppanggolin rarefaction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |

#### Optional arguments for ppanggolin rarefaction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-b, --beta` | float | 2.5 | beta is the strength of the smoothing using the graph topology during  partitioning. 0 will deactivate spatial smoothing. |
| `--depth` | int | 30 | Number of samplings at each sampling point |
| `--min` | int | 1 | Minimum number of genomes in a sample |
| `--max` | float | 100 | Maximum number of genomes in a sample (if above the number of provided genomes, the provided genomes will be the maximum) |
| `-ms, --max_degree_smoothing` | float | 10 | max. degree of the nodes to be included in the smoothing process. |
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `-fd, --free_dispersion` | bool | False | use if the dispersion around the centroid vector of each partition during must be free. It will be the same for all genomes by default. |
| `-ck, --chunk_size` | int | 500 | Size of the chunks when performing partitioning using chunks of genomes. Chunk partitioning will be used automatically if the number of genomes is above this number. |
| `-K, --nb_of_partitions` | int | -1 | Number of partitions to use. Must be at least 2. By default reuse K if it exists else compute it. |
| `--reestimate_K` | bool | False | Will recompute the number of partitions for each sample (between the values provided by --krange) (VERY intensive. Can take a long time.) |
| `-Kmm, --krange` | int | `[3, -1]` | Range of K values to test when detecting K automatically. Default between 3 and the K previously computed if there is one, or 20 if there are none. |
| `--soft_core` | float | 0.95 | Soft core threshold |
| `-se, --seed` | int | 42 | seed used to generate random numbers |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |

#### Common arguments for ppanggolin rarefaction

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


## Output

### `ppanggolin draw`

Draw figures representing the pangenome through different aspects

#### Required arguments for ppanggolin draw

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome.h5 file |

#### Optional arguments for ppanggolin draw

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `--tile_plot` | bool | False | draw the tile plot of the pangenome |
| `--soft_core` | float | 0.95 | Soft core threshold to use |
| `--ucurve` | bool | False | draw the U-curve of the pangenome |
| `--draw_spots` | bool | False | draw plots for spots of the pangenome |
| `--spots` | str | `all` | a comma-separated list of spots to draw (or 'all' to draw all spots, or 'synteny' to draw spots with different RGP syntenies). |
| `--nocloud` | bool | False | Do not draw the cloud genes in the tile plot |
| `--add_dendrogram` | bool | False | Include a dendrogram for genomes in the tile plot based on the presence/absence of gene families. |
| `--add_metadata` | bool | False | Display gene metadata as hover text for each cell in the tile plot. |
| `--metadata_sources` | str | — | Which source of metadata should be written in the tile plot. By default all metadata sources are included. |

#### Common arguments for ppanggolin draw

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin fasta`

Writes fasta files for different elements of the pangenome.

#### Required arguments for ppanggolin fasta

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Output directory where the file(s) will be written Required: Yes |

#### Contextually required arguments for ppanggolin fasta

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fasta` | Path | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed with gzip). One line per genome. |
| `--anno` | Path | — | A tab-separated file listing the genome names, and the gff/gbff filepath of its annotations (the files can be compressed with gzip). One line per genome. If this is provided, those annotations will be used. |

#### Output file for ppanggolin fasta

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--genes` | str | — | Write all nucleotide CDS sequences. Possible values are all, persistent, shell, cloud, rgp, softcore, core, module_X with X being a module id. |
| `--proteins` | str | — | Write representative amino acid sequences of genes. Possible values are all, persistent, shell, cloud, rgp, softcore, core, module_X with X being a module id. |
| `--prot_families` | str | — | Write representative amino acid sequences of gene families. Possible values are all, persistent, shell, cloud, rgp, softcore, core, module_X with X being a module id. |
| `--gene_families` | str | — | Write representative nucleotide sequences of gene families. Possible values are all, persistent, shell, cloud, rgp, softcore, core, module_X with X being a module id. |

#### Optional arguments for ppanggolin fasta

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--regions` | str | — | Write the RGP nucleotide sequences (requires --anno or --fasta used to compute the pangenome to be given) <br>Choices: `all`, `complete` |
| `--soft_core` | float | 0.95 | Soft core threshold to use if 'softcore' partition is chosen |
| `--compress` | bool | False | Compress the files in .gz |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `--cpu` | int | 1 | Number of available threads |
| `--tmpdir` | Path | `/tmp` | directory for storing temporary files |
| `--keep_tmp` | bool | False | Keeping temporary files (useful for debugging). |

#### Common arguments for ppanggolin fasta

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin info`

Prints information about a given pangenome graph file.

#### Required arguments for ppanggolin info

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | Path to the pangenome .h5 file Required: Yes |

#### Information Display Options (default: all) for ppanggolin info

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-a, --parameters` | bool | False | Display the parameters used or computed for each step of pangenome generation |
| `-c, --content` | bool | False | Display detailed information about the pangenome's content |
| `-s, --status` | bool | False | Display information about the statuses of different elements in the pangenome, indicating what has been computed or not |
| `-m, --metadata` | bool | False | Display a summary of the metadata saved in the pangenome |

#### Common arguments for ppanggolin info

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin metrics`

Compute several metrics on a given pangenome.

#### Required arguments for ppanggolin metrics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | Path to the pangenome .h5 file |

#### Input file for ppanggolin metrics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--genome_fluidity` | bool | False | Compute the pangenome genomic fluidity. |

#### Optional arguments for ppanggolin metrics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--no_print_info` | bool | False | Suppress printing the metrics result. Metrics are saved in the pangenome and viewable using 'ppanggolin info'. |
| `--recompute_metrics` | bool | False | Force re-computation of metrics if already computed. |

#### Common arguments for ppanggolin metrics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin write_genomes`

Writes 'flat' files that represent the genomes along with their associated pangenome elements.

#### Required arguments for ppanggolin write_genomes

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Output directory where the file(s) will be written Required: Yes |

#### Optional arguments for ppanggolin write_genomes

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--table` | bool | False | Generate a tsv file for each genome with pangenome annotations. |
| `--gff` | bool | False | Generate a gff file for each genome containing pangenome annotations. |
| `--proksee` | bool | False | Generate JSON map files for PROKSEE for each genome containing pangenome annotations to be used in proksee. |
| `--compress` | bool | False | Compress the files in .gz |
| `--genomes` | str | `all` | Specify the genomes for which to generate output. You can provide a list of genome names either directly in the command line separated by commas, or by referencing a file containing the list of genome names, with one name per line. |
| `--add_metadata` | bool | False | Include metadata information in the output files if any have been added to pangenome elements (see ppanggolin metadata command). |
| `--metadata_sources` | str | — | Which source of metadata should be written. By default all metadata sources are included. |
| `--metadata_sep` | str | `|` | The separator used to join multiple metadata values for elements with multiple metadata values from the same source. This character should not appear in metadata values. |
| `-c, --cpu` | int | 1 | Number of available cpus |

#### Contextually required arguments for ppanggolin write_genomes

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fasta` | Path | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed with gzip). One line per genome. |
| `--anno` | Path | — | A tab-separated file listing the genome names, and the gff/gbff filepath of its annotations (the files can be compressed with gzip). One line per genome. If this is provided, those annotations will be used. |

#### Common arguments for ppanggolin write_genomes

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin write_metadata`

Writes 'TSV' files that represent the metadata associated with elements of the pangenome.

#### Required arguments for ppanggolin write_metadata

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Output directory where the file(s) will be written Required: Yes |

#### Optional arguments for ppanggolin write_metadata

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--compress` | bool | False | Compress the files in .gz |
| `-e, --pangenome_elements` | set | `['contigs', 'families', 'genomes', 'modules', 'spots', 'RGPs', 'genes']` | Specify pangenome elements for which to write metadata. default is all element with metadata. <br>Choices: `contigs`, `families`, `genomes`, `modules`, `spots`, `RGPs`, `genes` |
| `-s, --metadata_sources` | str | — | Which source of metadata should be written. By default all metadata sources are included. |

#### Common arguments for ppanggolin write_metadata

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin write_pangenome`

Writes 'flat' files that represent the pangenome and its elements for use with other software.

#### Required arguments for ppanggolin write_pangenome

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Output directory where the file(s) will be written Required: Yes |

#### Optional arguments for ppanggolin write_pangenome

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--soft_core` | float | 0.95 | Soft core threshold to use |
| `--dup_margin` | float | 0.05 | Minimum ratio of genomes in which the family must have multiple genes for it to be considered 'duplicated' |
| `--gexf` | bool | False | Generate a detailed GEXF file with all genes and annotations for each family |
| `--light_gexf` | bool | False | Generate a simplified GEXF file with basic gene family information |
| `--gt` | bool | False | Generate a simplified graph-tool GT file with basic gene family and pangenome organism strains information |
| `--json` | bool | False | Export the graph in JSON file format |
| `--csv` | bool | False | Export gene presence/absence in CSV format (compatible with Roary). Uses partitions as alternative gene IDs if available. |
| `--Rtab` | bool | False | Export gene presence/absence as a tabular binary matrix |
| `--stats` | bool | False | Generate genome statistics ('genome_statistics.tsv') and duplication summary of persistent families ('mean_persistent_duplication.tsv') |
| `--partitions` | bool | False | Generate one file per partition listing the gene families it contains |
| `--families_tsv` | bool | False | Write a TSV file mapping genes to their corresponding gene families |
| `--regions` | bool | False | Write predicted RGPs (Regions of Genomic Plasticity) and metrics to 'plastic_regions.tsv' |
| `--regions_families` | bool | False | Write a TSV file mapping each RGP to its gene family content in 'rgp_families.tsv' |
| `--spots` | bool | False | Write spot summary and list all RGPs (Regions of Genomic Plasticity) per spot |
| `--borders` | bool | False | List all borders of each spot |
| `--modules` | bool | False | Write a TSV file listing functional modules and their associated gene families |
| `--spot_modules` | bool | False | Generate two files comparing module presence across spots |
| `--compress` | bool | False | Compress output files using gzip (.gz) |
| `-c, --cpu` | int | 1 | Number of available cpus |

#### Common arguments for ppanggolin write_pangenome

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


## Regions of Genomic Plasticity

### `ppanggolin module`

Predicts functional modules in your pangenome.

#### Required arguments for ppanggolin module

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |

#### Optional arguments for ppanggolin module

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--size` | int | 3 | Minimal number of gene family in a module |
| `--dup_margin` | float | 0.05 | minimum ratio of genomes in which the family must have multiple genes for it to be considered 'duplicated' |
| `-m, --min_presence` | int | 2 | Minimum number of times the module needs to be present in the pangenome to be reported. Increasing it will improve precision but lower sensitivity. |
| `-t, --transitive` | int | 4 | Size of the transitive closure used to build the graph. This indicates the number of non related genes allowed in-between two related genes. Increasing it will improve precision but lower sensitivity a little. |
| `-s, --jaccard` | float | 0.85 | minimum jaccard similarity used to filter edges between gene families. Increasing it will improve precision but lower sensitivity a lot. |
| `-c, --cpu` | int | 1 | Number of available cpus |

#### Common arguments for ppanggolin module

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin rgp`

Predicts Regions of Genomic Plasticity in the genomes of your pangenome.

#### Required arguments for ppanggolin rgp

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |

#### Optional arguments for ppanggolin rgp

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--persistent_penalty` | int | 3 | Penalty score to apply to persistent genes |
| `--variable_gain` | int | 1 | Gain score to apply to variable genes |
| `--min_score` | int | 4 | Minimal score wanted for considering a region as being a RGP |
| `--min_length` | int | 3000 | Minimum length (bp) of a region to be considered a RGP |
| `--dup_margin` | float | 0.05 | Minimum ratio of genomes where the family is present in which the family must have multiple genes for it to be considered 'duplicated' |

#### Common arguments for ppanggolin rgp

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin rgp_cluster`

Cluster RGPs based on their gene families.

#### Required arguments for ppanggolin rgp_cluster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file Required: Yes |

#### Optional arguments for ppanggolin rgp_cluster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--grr_cutoff` | float | 0.8 | Min gene repertoire relatedness metric used in the rgp clustering |
| `--grr_metric` | str | `incomplete_aware_grr` | The grr (Gene Repertoire Relatedness) is used to assess the similarity between two RGPs based on their gene families.There are three different modes for calculating the grr value: 'min_grr', 'max_grr' or  'incomplete_aware_grr'.'min_grr': Computes the number of gene families shared between the two RGPs and divides it by the smaller number of gene families among the two RGPs.'max_grr': Calculates the number of gene families shared between the two RGPs and divides it by the larger number of gene families among the two RGPs.'incomplete_aware_grr' (default): If at least one RGP is considered incomplete, which occurs when it is located at the border of a contig,the 'min_grr' mode is used. Otherwise, the 'max_grr' mode is applied. <br>Choices: `incomplete_aware_grr`, `min_grr`, `max_grr` |
| `--ignore_incomplete_rgp` | bool | False | Do not cluster RGPs located on a contig border which are likely incomplete. |
| `--no_identical_rgp_merging` | bool | False | Do not merge in one node identical RGP (i.e. having the same family content) before clustering. |
| `--basename` | str | `rgp_cluster` | basename for the output file |
| `-o, --output` | Path | `rgp_clustering` | Output directory |
| `--graph_formats` | str | `['gexf', 'graphml']` | Format of the output graph. <br>Choices: `gexf`, `graphml` |
| `--add_metadata` | bool | False | Include metadata information in the output files if any have been added to pangenome elements (see ppanggolin metadata command). |
| `--metadata_sources` | str | — | Which source of metadata should be written. By default all metadata sources are included. |
| `--metadata_sep` | str | `|` | The separator used to join multiple metadata values for elements with multiple metadata values from the same source. This character should not appear in metadata values. |

#### Common arguments for ppanggolin rgp_cluster

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin spot`

Predicts spots in your pangenome.

#### Required arguments for ppanggolin spot

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |

#### Optional arguments for ppanggolin spot

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `ppanggolin_output_<date>_<pid>` | Output directory |
| `--spot_graph` | bool | False | Writes a graph of pairs of blocks of single copy markers flanking RGPs, supposedly belonging to the same hotspot |
| `--overlapping_match` | int | 2 | The number of 'missing' persistent genes allowed when comparing flanking genes during hotspot computations |
| `--set_size` | int | 3 | Number of single copy markers to use as flanking genes for a RGP during hotspot computation |
| `--exact_match_size` | int | 1 | Number of perfectly matching flanking single copy markers required to associate RGPs during hotspot computation (Ex: If set to 1, two RGPs are in the same hotspot if both their 1st flanking genes are the same) |
| `--graph_formats` | str | `['gexf']` | Format of the output graph. <br>Choices: `gexf`, `graphml` |

#### Common arguments for ppanggolin spot

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


## Analysis using reference pangenomes

### `ppanggolin align`

Aligns a genome or a set of proteins to the pangenome gene families and predicts information from it.

#### Required arguments for ppanggolin align

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-S, --sequences` | Path | — | sequences (nucleotides or amino acids) to align on the pangenome gene families |
| `-p, --pangenome` | Path | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Output directory where the file(s) will be written Required: Yes |

#### Optional arguments for ppanggolin align

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--no_defrag` | bool | False | DO NOT Realign gene families to link fragments withtheir non-fragmented gene family. (default: False) |
| `--identity` | float | 0.5 | min identity percentage threshold |
| `--coverage` | float | 0.8 | min coverage percentage threshold |
| `--fast` | bool | False | Use representative sequences of gene families for input gene alignment. This option is faster but may be less sensitive. By default, all pangenome genes are used. |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `--getinfo` | bool | False | Use this option to extract info related to the best hit of each query, such as the RGP it is in, or the spots. |
| `--draw_related` | bool | False | Draw figures and provide graphs in a gexf format of the eventual spots associated to the input sequences |
| `--use_pseudo` | bool | False | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--tmpdir` | Path | `/tmp` | directory for storing temporary files |
| `--keep_tmp` | bool | False | Keeping temporary files (useful for debugging). |

#### Common arguments for ppanggolin align

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin context`

Local genomic context analysis.

#### Required arguments for ppanggolin context

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome.h5 file |
| `-o, --output` | Path | `ppanggolin_context_<date>_<pid>` | Output directory where the file(s) will be written |

#### Input file for ppanggolin context

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-S, --sequences` | Path | — | Fasta file with the sequences of interest |
| `-F, --family` | Path | — | List of family IDs of interest from the pangenome |

#### Optional arguments for ppanggolin context

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-t, --transitive` | int | 4 | Size of the transitive closure used to build the graph. This indicates the number of non related genes allowed in-between two related genes. Increasing it will improve precision but lower sensitivity a little. |
| `-w, --window_size` | int | 5 | Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. |
| `-s, --jaccard` | float | 0.85 | minimum jaccard similarity used to filter edges between gene families. Increasing it will improve precision but lower sensitivity a lot. |
| `--graph_format` | str | `graphml` | Format of the context graph. Can be gexf or graphml. <br>Choices: `gexf`, `graphml` |

#### Alignment arguments for ppanggolin context

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--no_defrag` | bool | False | DO NOT Realign gene families to link fragments withtheir non-fragmented gene family. |
| `--fast` | bool | False | Use representative sequences of gene families for input gene alignment. This option is recommended for faster processing but may be less sensitive. By default, all pangenome genes are used for alignment. |
| `--identity` | float | 0.8 | min identity percentage threshold |
| `--coverage` | float | 0.8 | min coverage percentage threshold |
| `--translation_table` | int | 11 | The translation table to use when the input sequences are nucleotide sequences. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |
| `--keep_tmp` | bool | False | Keeping temporary files (useful for debugging). |
| `-c, --cpu` | int | 1 | Number of available cpus |

#### Common arguments for ppanggolin context

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin msa`

Compute Multiple Sequence Alignments for pangenome gene families.

#### Required arguments for ppanggolin msa

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Output directory where the file(s) will be written Required: Yes |

#### Optional arguments. Indicating 'all' writes all elements. Writing a partition ('persistent', 'shell', 'cloud', 'core' or 'accessory') write the elements associated to said partition. for ppanggolin msa

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--soft_core` | float | 0.95 | Soft core threshold to use if 'softcore' partition is chosen |
| `--dup_margin` | float | 0.05 | minimum ratio of genomes in which the family must have multiple genes for it to be considered 'duplicated' |
| `--single_copy` | bool | False | Use report gene families that are considered 'single copy', for details see option --dup_margin |
| `--partition` | str | `core` | compute Multiple Sequence Alignment of the gene families in the given partition <br>Choices: `all`, `persistent`, `shell`, `cloud`, `core`, `accessory`, `softcore` |
| `--source` | str | `protein` | indicates whether to use protein or dna sequences to compute the msa <br>Choices: `dna`, `protein` |
| `--phylo` | bool | False | Writes a whole genome msa file for additional phylogenetic analysis |
| `--use_gene_id` | bool | False | Use gene identifiers rather than genome names for sequences in the family MSA (genome names are used by default) |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--tmpdir` | str | `/tmp` | directory for storing temporary files |

#### Common arguments for ppanggolin msa

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin projection`

Annotates external genomes with an existing pangenome.

#### Required arguments for ppanggolin projection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-p, --pangenome` | Path | — | The pangenome.h5 file |
| `--fasta` | Path | — | Specify a FASTA file containing the genomic sequences of the genome(s) you wish to annotate, or provide a tab-separated file listing genome names alongside their respective FASTA filepaths, with one line per genome. |
| `--anno` | Path | — | Specify an annotation file in GFF/GBFF format for the genome you wish to annotate. Alternatively, you can provide a tab-separated file listing genome names alongside their respective annotation filepaths, with one line per genome. If both an annotation file and a FASTA file are provided, the annotation file will take precedence. |

#### Single Genome Arguments for ppanggolin projection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-n, --genome_name` | str | `input_genome` | Specify the name of the genome whose genome you want to annotate when providing a single FASTA or annotation file. |
| `--circular_contigs` | tuple | — | Specify the contigs of the input genome that should be treated as circular when providing a single FASTA or annotation file. |

#### Optional arguments for ppanggolin projection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `ppanggolin_projection_<date>_<pid>` | Output directory |
| `--translation_table` | int | 11 | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `--no_defrag` | bool | False | DO NOT Realign gene families to link fragments with their non-fragmented gene family. (default: False) |
| `--fast` | bool | False | Use representative sequences of gene families for input gene alignment. This option is faster but may be less sensitive. By default, all pangenome genes are used. |
| `--identity` | float | 0.8 | min identity percentage threshold |
| `--coverage` | float | 0.8 | min coverage percentage threshold |
| `--use_pseudo` | bool | False | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `--dup_margin` | float | 0.05 | minimum ratio of genomes in which the family must have multiple genes for it to be considered 'duplicated'. This metric is used to compute completeness and duplication of the input genomes |
| `--soft_core` | float | 0.95 | Soft core threshold used when generating general statistics on the projected genome. This threshold does not influence PPanGGOLiN's partitioning. The value determines the minimum fraction of genomes that must possess a gene family for it to be considered part of the soft core. |
| `--spot_graph` | bool | False | Write the spot graph to a file, with pairs of blocks of single copy markers flanking RGPs as nodes. This graph can be used to visualize nodes that have RGPs from the input genome. |
| `--graph_formats` | str | `['gexf']` | Format of the output graph. <br>Choices: `gexf`, `graphml` |
| `--gff` | bool | False | Generate GFF files with projected pangenome annotations for each input genome. |
| `--proksee` | bool | False | Generate JSON map files for PROKSEE with projected pangenome annotations for each input genome. |
| `--table` | bool | False | Generate a tsv file for each input genome with pangenome annotations. |
| `--compress` | bool | False | Compress the files in .gz |
| `--add_sequences` | bool | False | Include input genome DNA sequences in GFF and Proksee output. |
| `-c, --cpu` | int | 1 | Number of available cpus |
| `--tmpdir` | Path | `/tmp` | directory for storing temporary files |
| `--keep_tmp` | bool | False | Keeping temporary files (useful for debugging). |
| `--add_metadata` | bool | False | Include metadata information in the output files if any have been added to pangenome elements (see ppanggolin metadata command). |
| `--metadata_sources` | str | — | Which source of metadata should be written. By default all metadata sources are included. |
| `--metadata_sep` | str | `|` | The separator used to join multiple metadata values for elements with multiple metadata values from the same source. This character should not appear in metadata values. |

#### Common arguments for ppanggolin projection

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |


## Utility command

### `ppanggolin utils`

Helper side commands.

#### Required arguments for ppanggolin utils

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--default_config` | str | — | Generate a config file with default values for the given subcommand. <br>Choices: `annotate`, `cluster`, `graph`, `partition`, `rarefaction`, `workflow`, `panrgp`, `panmodule`, `all`, `draw`, `write_pangenome`, `write_genomes`, `write_metadata`, `fasta`, `msa`, `metrics`, `align`, `info`, `rgp`, `spot`, `module`, `context`, `projection`, `rgp_cluster`, `metadata`, `utils` |

#### Config arguments for ppanggolin utils

| Parameter | Type | Default | Description |
|---|---|---|---|
| `-o, --output` | Path | `default_config.yaml` | name and path of the config file with default parameters written in yaml. |

#### Common arguments for ppanggolin utils

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--verbose` | int | 1 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) <br>Choices: `0`, `1`, `2` |
| `--log` | str | `stdout` | log output file |
| `-d, --disable_prog_bar` | bool | False | disables the progress bars |
| `-f, --force` | bool | False | Force writing in output directory and in pangenome output file. |
| `--config` | Path | — | Specify command arguments through a YAML configuration file. |

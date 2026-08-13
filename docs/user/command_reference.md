# Command-line reference

Complete reference of the PPanGGOLiN command-line interface.

This page is automatically generated from the PPanGGOLiN command-line parser.

## Basic

### `ppanggolin all`

Easy workflow to run all possible analysis

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `--fasta` | Path | — | No | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per genome. This option can be used alone. |
| `--anno` | Path | — | No | — | A tab-separated file listing the genome names, and the gff filepath of its annotations (the gffs can be compressed). One line per genome. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used. |
| `--clusters` | Path | — | No | — | a tab-separated file listing the cluster names, the gene IDs, and optionally whether they are a fragment or not. |
| `-o, --output` | Path | ppanggolin_output_DATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `--basename` | str | `pangenome` | No | — | basename for the output file |
| `--rarefaction` | 0 | False | No | — | Use to compute the rarefaction curves (WARNING: can be time consuming) |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--kingdom` | lower | `bacteria` | No | bacteria, archaea | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. |
| `--mode` | str | `1` | No | 0, 1, 2, 3 | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) |
| `--coverage` | restricted_float | 0.8 | No | — | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--identity` | restricted_float | 0.8 | No | — | Minimal identity percent for two proteins to be in the same cluster |
| `--infer_singletons` | 0 | False | No | — | Use this option together with --clusters. If a gene is not present in the provided clustering result file, it will be assigned to its own unique cluster as a singleton. |
| `--use_pseudo` | 0 | False | No | — | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-K, --nb_of_partitions` | int | -1 | No | — | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `--no_defrag` | 0 | False | No | — | DO NOT Realign gene families to link fragments with their non-fragmented gene family. |
| `--no_flat_files` | 0 | False | No | — | Generate only the HDF5 pangenome file. |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin panmodule`

Easy workflow to run a pangenome analysis with module prediction

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `--fasta` | Path | — | No | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per genome. This option can be used alone. |
| `--anno` | Path | — | No | — | A tab-separated file listing the genome names, and the gff filepath of its annotations (the gffs can be compressed). One line per genome. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used. |
| `--clusters` | Path | — | No | — | a tab-separated file listing the cluster names, the gene IDs, and optionally whether they are a fragment or not. |
| `-o, --output` | Path | ppanggolin_output_DATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `--basename` | str | `pangenome` | No | — | basename for the output file |
| `--rarefaction` | 0 | False | No | — | Use to compute the rarefaction curves (WARNING: can be time consuming) |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--kingdom` | lower | `bacteria` | No | bacteria, archaea | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. |
| `--mode` | str | `1` | No | 0, 1, 2, 3 | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) |
| `--coverage` | restricted_float | 0.8 | No | — | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--identity` | restricted_float | 0.8 | No | — | Minimal identity percent for two proteins to be in the same cluster |
| `--infer_singletons` | 0 | False | No | — | Use this option together with --clusters. If a gene is not present in the provided clustering result file, it will be assigned to its own unique cluster as a singleton. |
| `--use_pseudo` | 0 | False | No | — | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-K, --nb_of_partitions` | int | -1 | No | — | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `--no_defrag` | 0 | False | No | — | DO NOT Realign gene families to link fragments with their non-fragmented gene family. |
| `--no_flat_files` | 0 | False | No | — | Generate only the HDF5 pangenome file. |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin panrgp`

Easy workflow to run a pangenome analysis with genomic islands and spots of insertion detection

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `--fasta` | Path | — | No | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per genome. This option can be used alone. |
| `--anno` | Path | — | No | — | A tab-separated file listing the genome names, and the gff filepath of its annotations (the gffs can be compressed). One line per genome. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used. |
| `--clusters` | Path | — | No | — | a tab-separated file listing the cluster names, the gene IDs, and optionally whether they are a fragment or not. |
| `-o, --output` | Path | ppanggolin_output_DATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `--basename` | str | `pangenome` | No | — | basename for the output file |
| `--rarefaction` | 0 | False | No | — | Use to compute the rarefaction curves (WARNING: can be time consuming) |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--kingdom` | lower | `bacteria` | No | bacteria, archaea | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. |
| `--mode` | str | `1` | No | 0, 1, 2, 3 | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) |
| `--coverage` | restricted_float | 0.8 | No | — | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--identity` | restricted_float | 0.8 | No | — | Minimal identity percent for two proteins to be in the same cluster |
| `--infer_singletons` | 0 | False | No | — | Use this option together with --clusters. If a gene is not present in the provided clustering result file, it will be assigned to its own unique cluster as a singleton. |
| `--use_pseudo` | 0 | False | No | — | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-K, --nb_of_partitions` | int | -1 | No | — | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `--no_defrag` | 0 | False | No | — | DO NOT Realign gene families to link fragments with their non-fragmented gene family. |
| `--no_flat_files` | 0 | False | No | — | Generate only the HDF5 pangenome file. |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin workflow`

Easy workflow to run a pangenome analysis in one go

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `--fasta` | Path | — | No | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per genome. This option can be used alone. |
| `--anno` | Path | — | No | — | A tab-separated file listing the genome names, and the gff filepath of its annotations (the gffs can be compressed). One line per genome. This option can be used alone IF the fasta sequences are in the gff files, otherwise --fasta needs to be used. |
| `--clusters` | Path | — | No | — | a tab-separated file listing the cluster names, the gene IDs, and optionally whether they are a fragment or not. |
| `-o, --output` | Path | ppanggolin_output_DATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `--basename` | str | `pangenome` | No | — | basename for the output file |
| `--rarefaction` | 0 | False | No | — | Use to compute the rarefaction curves (WARNING: can be time consuming) |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--kingdom` | lower | `bacteria` | No | bacteria, archaea | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. |
| `--mode` | str | `1` | No | 0, 1, 2, 3 | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) |
| `--coverage` | restricted_float | 0.8 | No | — | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--identity` | restricted_float | 0.8 | No | — | Minimal identity percent for two proteins to be in the same cluster |
| `--infer_singletons` | 0 | False | No | — | Use this option together with --clusters. If a gene is not present in the provided clustering result file, it will be assigned to its own unique cluster as a singleton. |
| `--use_pseudo` | 0 | False | No | — | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-K, --nb_of_partitions` | int | -1 | No | — | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `--no_defrag` | 0 | False | No | — | DO NOT Realign gene families to link fragments with their non-fragmented gene family. |
| `--no_flat_files` | 0 | False | No | — | Generate only the HDF5 pangenome file. |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


## Expert

### `ppanggolin annotate`

Annotate genomes

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `--fasta` | Path | — | No | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed with gzip). One line per genome. |
| `--anno` | Path | — | No | — | A tab-separated file listing the genome names, and the gff/gbff filepath of its annotations (the files can be compressed with gzip). One line per genome. If this is provided, those annotations will be used. |
| `-o, --output` | Path | ppanggolin_output_DATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `--allow_overlap` | 0 | False | No | — | Use to not remove genes overlapping with RNA features. |
| `--norna` | 0 | False | No | — | Use to avoid annotating RNA features. |
| `--kingdom` | lower | `bacteria` | No | bacteria, archaea | Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation. |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified and using annotation files (--anno), the translation table information found in the annotation files will be used. Otherwise, the default genetic code 11 will be used. |
| `--basename` | str | `pangenome` | No | — | basename for the output file |
| `--use_pseudo` | 0 | False | No | — | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-p, --prodigal_procedure` | lower | — | No | single, meta | Allow to force the prodigal procedure. If nothing given, PPanGGOLiN will decide in function of contig length |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin cluster`

Cluster genes into gene families

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `--identity` | restricted_float | 0.8 | No | — | Minimal identity percent for two proteins to be in the same cluster |
| `--coverage` | restricted_float | 0.8 | No | — | Minimal coverage of the alignment for two proteins to be in the same cluster |
| `--mode` | str | `1` | No | 0, 1, 2, 3 | the cluster mode of MMseqs2. 0: Setcover, 1: single linkage (or connected component), 2: CD-HIT-like, 3: CD-HIT-like (lowmem) |
| `--no_defrag` | 0 | False | No | — | DO NOT Use the defragmentation strategy to link potential fragments with their original gene family. |
| `--clusters` | Path | — | No | — | A tab-separated list containing the result of a clustering. One line per gene. First column is cluster ID, and second is gene ID |
| `--infer_singletons` | 0 | False | No | — | When reading a clustering result with --clusters, if a gene is not in the provided file it will be placed in a cluster where the gene is the only member. |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--tmpdir` | Path | /tmp | No | — | directory for storing temporary files |
| `--keep_tmp` | 0 | False | No | — | Keeping temporary files (useful for debugging). |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin graph`

Create the pangenome graph

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-r, --remove_high_copy_number` | int | 0 | No | — | Positive Number: Remove families having a number of copy of gene in a single genome above or equal to this threshold in at least one genome (0 or negative values are ignored). |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin metadata`

Add metadata to elements in yout pangenome

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-m, --metadata` | Path | — | No | — | Metadata in TSV file. See our github for more detail about format |
| `-s, --source` | str | — | No | — | Name of the metadata source |
| `-a, --assign` | str | — | No | families, genomes, contigs, genes, RGPs, spots, modules | Select to which pangenome element metadata will be assigned |
| `--omit` | 0 | False | No | — | Allow to pass if a key in metadata is not find in pangenome |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin partition`

Partition the pangenome graph

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome.h5 file |
| `-b, --beta` | float | 2.5 | No | — | beta is the strength of the smoothing using the graph topology during partitioning. 0 will deactivate spatial smoothing. |
| `-ms, --max_degree_smoothing` | float | 10 | No | — | max. degree of the nodes to be included in the smoothing process. |
| `-o, --output` | Path | ppanggolin_outputDATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `-fd, --free_dispersion` | 0 | False | No | — | use if the dispersion around the centroid vector of each partition during must be free. It will be the same for all genomes by default. |
| `-ck, --chunk_size` | int | 500 | No | — | Size of the chunks when performing partitioning using chunks of genomes. Chunk partitioning will be used automatically if the number of genomes is above this number. |
| `-K, --nb_of_partitions` | int | -1 | No | — | Number of partitions to use. Must be at least 2. If under 2, it will be detected automatically. |
| `-Kmm, --krange` | 2 | [3, 20] | No | — | Range of K values to test when detecting K automatically. |
| `-im, --ICL_margin` | float | 0.05 | No | — | K is detected automatically by maximizing ICL. However at some point the ICL reaches a plateau. Therefore we are looking for the minimal value of K without significant gain from the larger values of K measured by ICL. For that we take the lowest K that is found within a given 'margin' of the maximal ICL value. Basically, change this option only if you truly understand it, otherwise just leave it be. |
| `--draw_ICL` | 0 | False | No | — | Use if you want to draw the ICL curve for all the tested K values. Will not be done if K is given. |
| `--keep_tmp_files` | 0 | False | No | — | Use if you want to keep the temporary NEM files |
| `-se, --seed` | int | 42 | No | — | seed used to generate random numbers |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin rarefaction`

Compute the rarefaction curve of the pangenome

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-b, --beta` | float | 2.5 | No | — | beta is the strength of the smoothing using the graph topology during  partitioning. 0 will deactivate spatial smoothing. |
| `--depth` | int | 30 | No | — | Number of samplings at each sampling point |
| `--min` | int | 1 | No | — | Minimum number of genomes in a sample |
| `--max` | float | 100 | No | — | Maximum number of genomes in a sample (if above the number of provided genomes, the provided genomes will be the maximum) |
| `-ms, --max_degree_smoothing` | float | 10 | No | — | max. degree of the nodes to be included in the smoothing process. |
| `-o, --output` | Path | ppanggolin_outputDATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `-fd, --free_dispersion` | 0 | False | No | — | use if the dispersion around the centroid vector of each partition during must be free. It will be the same for all genomes by default. |
| `-ck, --chunk_size` | int | 500 | No | — | Size of the chunks when performing partitioning using chunks of genomes. Chunk partitioning will be used automatically if the number of genomes is above this number. |
| `-K, --nb_of_partitions` | int | -1 | No | — | Number of partitions to use. Must be at least 2. By default reuse K if it exists else compute it. |
| `--reestimate_K` | 0 | False | No | — | Will recompute the number of partitions for each sample (between the values provided by --krange) (VERY intensive. Can take a long time.) |
| `-Kmm, --krange` | 2 | [3, -1] | No | — | Range of K values to test when detecting K automatically. Default between 3 and the K previously computed if there is one, or 20 if there are none. |
| `--soft_core` | float | 0.95 | No | — | Soft core threshold |
| `-se, --seed` | int | 42 | No | — | seed used to generate random numbers |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


## Output

### `ppanggolin draw`

Draw figures representing the pangenome through different aspects

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome.h5 file |
| `-o, --output` | Path | ppanggolin_output_DATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `--tile_plot` | 0 | False | No | — | draw the tile plot of the pangenome |
| `--soft_core` | restricted_float | 0.95 | No | — | Soft core threshold to use |
| `--ucurve` | 0 | False | No | — | draw the U-curve of the pangenome |
| `--draw_spots` | 0 | False | No | — | draw plots for spots of the pangenome |
| `--spots` | list | `all` | No | — | a comma-separated list of spots to draw (or 'all' to draw all spots, or 'synteny' to draw spots with different RGP syntenies). |
| `--nocloud` | 0 | False | No | — | Do not draw the cloud genes in the tile plot |
| `--add_dendrogram` | 0 | False | No | — | Include a dendrogram for genomes in the tile plot based on the presence/absence of gene families. |
| `--add_metadata` | 0 | False | No | — | Display gene metadata as hover text for each cell in the tile plot. |
| `--metadata_sources` | list | — | No | — | Which source of metadata should be written in the tile plot. By default all metadata sources are included. |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin fasta`

Writes fasta files for different elements of the pangenome.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Yes | — | Output directory where the file(s) will be written |
| `--fasta` | Path | — | No | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed with gzip). One line per genome. |
| `--anno` | Path | — | No | — | A tab-separated file listing the genome names, and the gff/gbff filepath of its annotations (the files can be compressed with gzip). One line per genome. If this is provided, those annotations will be used. |
| `--genes` | filter_values | — | No | — | Write all nucleotide CDS sequences. Possible values are all, persistent, shell, cloud, rgp, softcore, core, module_X with X being a module id. |
| `--proteins` | filter_values | — | No | — | Write representative amino acid sequences of genes. Possible values are all, persistent, shell, cloud, rgp, softcore, core, module_X with X being a module id. |
| `--prot_families` | filter_values | — | No | — | Write representative amino acid sequences of gene families. Possible values are all, persistent, shell, cloud, rgp, softcore, core, module_X with X being a module id. |
| `--gene_families` | filter_values | — | No | — | Write representative nucleotide sequences of gene families. Possible values are all, persistent, shell, cloud, rgp, softcore, core, module_X with X being a module id. |
| `--regions` | str | — | No | all, complete | Write the RGP nucleotide sequences (requires --anno or --fasta used to compute the pangenome to be given) |
| `--soft_core` | restricted_float | 0.95 | No | — | Soft core threshold to use if 'softcore' partition is chosen |
| `--compress` | 0 | False | No | — | Compress the files in .gz |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `--cpu` | int | 1 | No | — | Number of available threads |
| `--tmpdir` | Path | /tmp | No | — | directory for storing temporary files |
| `--keep_tmp` | 0 | False | No | — | Keeping temporary files (useful for debugging). |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin metrics`

Compute several metrics on a given pangenome.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | Path to the pangenome .h5 file |
| `--genome_fluidity` | 0 | False | No | — | Compute the pangenome genomic fluidity. |
| `--no_print_info` | 0 | False | No | — | Suppress printing the metrics result. Metrics are saved in the pangenome and viewable using 'ppanggolin info'. |
| `--recompute_metrics` | 0 | False | No | — | Force re-computation of metrics if already computed. |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin write_genomes`

Writes 'flat' files that represent the genomes along with their associated pangenome elements.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Yes | — | Output directory where the file(s) will be written |
| `--table` | 0 | False | No | — | Generate a tsv file for each genome with pangenome annotations. |
| `--gff` | 0 | False | No | — | Generate a gff file for each genome containing pangenome annotations. |
| `--proksee` | 0 | False | No | — | Generate JSON map files for PROKSEE for each genome containing pangenome annotations to be used in proksee. |
| `--compress` | 0 | False | No | — | Compress the files in .gz |
| `--genomes` | str | `all` | No | — | Specify the genomes for which to generate output. You can provide a list of genome names either directly in the command line separated by commas, or by referencing a file containing the list of genome names, with one name per line. |
| `--add_metadata` | 0 | False | No | — | Include metadata information in the output files if any have been added to pangenome elements (see ppanggolin metadata command). |
| `--metadata_sources` | list | — | No | — | Which source of metadata should be written. By default all metadata sources are included. |
| `--metadata_sep` | str | `|` | No | — | The separator used to join multiple metadata values for elements with multiple metadata values from the same source. This character should not appear in metadata values. |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--fasta` | Path | — | No | — | A tab-separated file listing the genome names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed with gzip). One line per genome. |
| `--anno` | Path | — | No | — | A tab-separated file listing the genome names, and the gff/gbff filepath of its annotations (the files can be compressed with gzip). One line per genome. If this is provided, those annotations will be used. |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin write_metadata`

Writes 'TSV' files that represent the metadata associated with elements of the pangenome.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Yes | — | Output directory where the file(s) will be written |
| `--compress` | 0 | False | No | — | Compress the files in .gz |
| `-e, --pangenome_elements` | list | ['families', 'genes', 'spots', 'genomes', 'RGPs', 'modules', 'contigs'] | No | families, genes, spots, genomes, RGPs, modules, contigs | Specify pangenome elements for which to write metadata. default is all element with metadata. |
| `-s, --metadata_sources` | list | — | No | — | Which source of metadata should be written. By default all metadata sources are included. |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin write_pangenome`

Writes 'flat' files that represent the pangenome and its elements for use with other software.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Yes | — | Output directory where the file(s) will be written |
| `--soft_core` | restricted_float | 0.95 | No | — | Soft core threshold to use |
| `--dup_margin` | restricted_float | 0.05 | No | — | Minimum ratio of genomes in which the family must have multiple genes for it to be considered 'duplicated' |
| `--gexf` | 0 | False | No | — | Generate a detailed GEXF file with all genes and annotations for each family |
| `--light_gexf` | 0 | False | No | — | Generate a simplified GEXF file with basic gene family information |
| `--gt` | 0 | False | No | — | Generate a simplified graph-tool GT file with basic gene family and pangenome organism strains information |
| `--json` | 0 | False | No | — | Export the graph in JSON file format |
| `--csv` | 0 | False | No | — | Export gene presence/absence in CSV format (compatible with Roary). Uses partitions as alternative gene IDs if available. |
| `--Rtab` | 0 | False | No | — | Export gene presence/absence as a tabular binary matrix |
| `--stats` | 0 | False | No | — | Generate genome statistics ('genome_statistics.tsv') and duplication summary of persistent families ('mean_persistent_duplication.tsv') |
| `--partitions` | 0 | False | No | — | Generate one file per partition listing the gene families it contains |
| `--families_tsv` | 0 | False | No | — | Write a TSV file mapping genes to their corresponding gene families |
| `--regions` | 0 | False | No | — | Write predicted RGPs (Regions of Genomic Plasticity) and metrics to 'plastic_regions.tsv' |
| `--regions_families` | 0 | False | No | — | Write a TSV file mapping each RGP to its gene family content in 'rgp_families.tsv' |
| `--spots` | 0 | False | No | — | Write spot summary and list all RGPs (Regions of Genomic Plasticity) per spot |
| `--borders` | 0 | False | No | — | List all borders of each spot |
| `--modules` | 0 | False | No | — | Write a TSV file listing functional modules and their associated gene families |
| `--spot_modules` | 0 | False | No | — | Generate two files comparing module presence across spots |
| `--compress` | 0 | False | No | — | Compress output files using gzip (.gz) |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


## Regions of Genomic Plasticity

### `ppanggolin module`

Predicts functional modules in your pangenome.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `--size` | int | 3 | No | — | Minimal number of gene family in a module |
| `--dup_margin` | restricted_float | 0.05 | No | — | minimum ratio of genomes in which the family must have multiple genes for it to be considered 'duplicated' |
| `-m, --min_presence` | int | 2 | No | — | Minimum number of times the module needs to be present in the pangenome to be reported. Increasing it will improve precision but lower sensitivity. |
| `-t, --transitive` | int | 4 | No | — | Size of the transitive closure used to build the graph. This indicates the number of non related genes allowed in-between two related genes. Increasing it will improve precision but lower sensitivity a little. |
| `-s, --jaccard` | restricted_float | 0.85 | No | — | minimum jaccard similarity used to filter edges between gene families. Increasing it will improve precision but lower sensitivity a lot. |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin rgp`

Predicts Regions of Genomic Plasticity in the genomes of your pangenome.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `--persistent_penalty` | int | 3 | No | — | Penalty score to apply to persistent genes |
| `--variable_gain` | int | 1 | No | — | Gain score to apply to variable genes |
| `--min_score` | int | 4 | No | — | Minimal score wanted for considering a region as being a RGP |
| `--min_length` | int | 3000 | No | — | Minimum length (bp) of a region to be considered a RGP |
| `--dup_margin` | restricted_float | 0.05 | No | — | Minimum ratio of genomes where the family is present in which the family must have multiple genes for it to be considered 'duplicated' |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin rgp_cluster`

Cluster RGPs based on their gene families.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | Yes | — | The pangenome .h5 file |
| `--grr_cutoff` | restricted_float | 0.8 | No | — | Min gene repertoire relatedness metric used in the rgp clustering |
| `--grr_metric` | str | `incomplete_aware_grr` | No | incomplete_aware_grr, min_grr, max_grr | The grr (Gene Repertoire Relatedness) is used to assess the similarity between two RGPs based on their gene families.There are three different modes for calculating the grr value: 'min_grr', 'max_grr' or  'incomplete_aware_grr'.'min_grr': Computes the number of gene families shared between the two RGPs and divides it by the smaller number of gene families among the two RGPs.'max_grr': Calculates the number of gene families shared between the two RGPs and divides it by the larger number of gene families among the two RGPs.'incomplete_aware_grr' (default): If at least one RGP is considered incomplete, which occurs when it is located at the border of a contig,the 'min_grr' mode is used. Otherwise, the 'max_grr' mode is applied. |
| `--ignore_incomplete_rgp` | 0 | False | No | — | Do not cluster RGPs located on a contig border which are likely incomplete. |
| `--no_identical_rgp_merging` | 0 | False | No | — | Do not merge in one node identical RGP (i.e. having the same family content) before clustering. |
| `--basename` | str | `rgp_cluster` | No | — | basename for the output file |
| `-o, --output` | Path | `rgp_clustering` | No | — | Output directory |
| `--graph_formats` | list | ['gexf', 'graphml'] | No | gexf, graphml | Format of the output graph. |
| `--add_metadata` | 0 | False | No | — | Include metadata information in the output files if any have been added to pangenome elements (see ppanggolin metadata command). |
| `--metadata_sources` | list | — | No | — | Which source of metadata should be written. By default all metadata sources are included. |
| `--metadata_sep` | str | `|` | No | — | The separator used to join multiple metadata values for elements with multiple metadata values from the same source. This character should not appear in metadata values. |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin spot`

Predicts spots in your pangenome.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-o, --output` | Path | ppanggolin_outputDATE2026-08-13_HOUR12.55.19_PID729243 | No | — | Output directory |
| `--spot_graph` | 0 | False | No | — | Writes a graph of pairs of blocks of single copy markers flanking RGPs, supposedly belonging to the same hotspot |
| `--overlapping_match` | int | 2 | No | — | The number of 'missing' persistent genes allowed when comparing flanking genes during hotspot computations |
| `--set_size` | int | 3 | No | — | Number of single copy markers to use as flanking genes for a RGP during hotspot computation |
| `--exact_match_size` | int | 1 | No | — | Number of perfectly matching flanking single copy markers required to associate RGPs during hotspot computation (Ex: If set to 1, two RGPs are in the same hotspot if both their 1st flanking genes are the same) |
| `--graph_formats` | list | ['gexf'] | No | gexf, graphml | Format of the output graph. |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


## Analysis using reference pangenomes

### `ppanggolin align`

Aligns a genome or a set of proteins to the pangenome gene families and predicts information from it.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-S, --sequences` | Path | — | No | — | sequences (nucleotides or amino acids) to align on the pangenome gene families |
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Yes | — | Output directory where the file(s) will be written |
| `--no_defrag` | 0 | False | No | — | DO NOT Realign gene families to link fragments withtheir non-fragmented gene family. (default: False) |
| `--identity` | float | 0.5 | No | — | min identity percentage threshold |
| `--coverage` | float | 0.8 | No | — | min coverage percentage threshold |
| `--fast` | 0 | False | No | — | Use representative sequences of gene families for input gene alignment. This option is faster but may be less sensitive. By default, all pangenome genes are used. |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `--getinfo` | 0 | False | No | — | Use this option to extract info related to the best hit of each query, such as the RGP it is in, or the spots. |
| `--draw_related` | 0 | False | No | — | Draw figures and provide graphs in a gexf format of the eventual spots associated to the input sequences |
| `--use_pseudo` | 0 | False | No | — | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--tmpdir` | Path | /tmp | No | — | directory for storing temporary files |
| `--keep_tmp` | 0 | False | No | — | Keeping temporary files (useful for debugging). |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin context`

Local genomic context analysis.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome.h5 file |
| `-o, --output` | Path | `ppanggolin_context_DATE2026-08-13_HOUR12.55.19_PID729243` | No | — | Output directory where the file(s) will be written |
| `-S, --sequences` | Path | — | No | — | Fasta file with the sequences of interest |
| `-F, --family` | Path | — | No | — | List of family IDs of interest from the pangenome |
| `-t, --transitive` | int | 4 | No | — | Size of the transitive closure used to build the graph. This indicates the number of non related genes allowed in-between two related genes. Increasing it will improve precision but lower sensitivity a little. |
| `-w, --window_size` | int | 5 | No | — | Number of neighboring genes that are considered on each side of a gene of interest when searching for conserved genomic contexts. |
| `-s, --jaccard` | restricted_float | 0.85 | No | — | minimum jaccard similarity used to filter edges between gene families. Increasing it will improve precision but lower sensitivity a lot. |
| `--graph_format` | str | `graphml` | No | gexf, graphml | Format of the context graph. Can be gexf or graphml. |
| `--no_defrag` | 0 | False | No | — | DO NOT Realign gene families to link fragments withtheir non-fragmented gene family. |
| `--fast` | 0 | False | No | — | Use representative sequences of gene families for input gene alignment. This option is recommended for faster processing but may be less sensitive. By default, all pangenome genes are used for alignment. |
| `--identity` | float | 0.8 | No | — | min identity percentage threshold |
| `--coverage` | float | 0.8 | No | — | min coverage percentage threshold |
| `--translation_table` | int | 11 | No | — | The translation table to use when the input sequences are nucleotide sequences. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--keep_tmp` | 0 | False | No | — | Keeping temporary files (useful for debugging). |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin msa`

Compute Multiple Sequence Alignments for pangenome gene families.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome .h5 file |
| `-o, --output` | Path | — | Yes | — | Output directory where the file(s) will be written |
| `--soft_core` | restricted_float | 0.95 | No | — | Soft core threshold to use if 'softcore' partition is chosen |
| `--dup_margin` | restricted_float | 0.05 | No | — | minimum ratio of genomes in which the family must have multiple genes for it to be considered 'duplicated' |
| `--single_copy` | 0 | False | No | — | Use report gene families that are considered 'single copy', for details see option --dup_margin |
| `--partition` | str | `core` | No | all, persistent, shell, cloud, core, accessory, softcore | compute Multiple Sequence Alignment of the gene families in the given partition |
| `--source` | str | `protein` | No | dna, protein | indicates whether to use protein or dna sequences to compute the msa |
| `--phylo` | 0 | False | No | — | Writes a whole genome msa file for additional phylogenetic analysis |
| `--use_gene_id` | 0 | False | No | — | Use gene identifiers rather than genome names for sequences in the family MSA (genome names are used by default) |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--tmpdir` | str | /tmp | No | — | directory for storing temporary files |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |


### `ppanggolin projection`

Annotates external genomes with an existing pangenome.

#### Parameters

| Parameter | Type | Default | Required | Choices | Description |
|---|---|---|---|---|---|
| `-p, --pangenome` | Path | — | No | — | The pangenome.h5 file |
| `--fasta` | Path | — | No | — | Specify a FASTA file containing the genomic sequences of the genome(s) you wish to annotate, or provide a tab-separated file listing genome names alongside their respective FASTA filepaths, with one line per genome. |
| `--anno` | Path | — | No | — | Specify an annotation file in GFF/GBFF format for the genome you wish to annotate. Alternatively, you can provide a tab-separated file listing genome names alongside their respective annotation filepaths, with one line per genome. If both an annotation file and a FASTA file are provided, the annotation file will take precedence. |
| `-n, --genome_name` | str | `input_genome` | No | — | Specify the name of the genome whose genome you want to annotate when providing a single FASTA or annotation file. |
| `--circular_contigs` | list | — | No | — | Specify the contigs of the input genome that should be treated as circular when providing a single FASTA or annotation file. |
| `-o, --output` | Path | `ppanggolin_projection_DATE2026-08-13_HOUR12.55.19_PID729243` | No | — | Output directory |
| `--translation_table` | int | 11 | No | — | Translation table (genetic code) to use. If not specified, the translation table used when building the pangenome will be used. This can be accessed using 'ppanggolin info'. |
| `--no_defrag` | 0 | False | No | — | DO NOT Realign gene families to link fragments with their non-fragmented gene family. (default: False) |
| `--fast` | 0 | False | No | — | Use representative sequences of gene families for input gene alignment. This option is faster but may be less sensitive. By default, all pangenome genes are used. |
| `--identity` | restricted_float | 0.8 | No | — | min identity percentage threshold |
| `--coverage` | restricted_float | 0.8 | No | — | min coverage percentage threshold |
| `--use_pseudo` | 0 | False | No | — | In the context of provided annotation, use this option to read pseudogenes. (Default behavior is to ignore them) |
| `--dup_margin` | restricted_float | 0.05 | No | — | minimum ratio of genomes in which the family must have multiple genes for it to be considered 'duplicated'. This metric is used to compute completeness and duplication of the input genomes |
| `--soft_core` | restricted_float | 0.95 | No | — | Soft core threshold used when generating general statistics on the projected genome. This threshold does not influence PPanGGOLiN's partitioning. The value determines the minimum fraction of genomes that must possess a gene family for it to be considered part of the soft core. |
| `--spot_graph` | 0 | False | No | — | Write the spot graph to a file, with pairs of blocks of single copy markers flanking RGPs as nodes. This graph can be used to visualize nodes that have RGPs from the input genome. |
| `--graph_formats` | list | ['gexf'] | No | gexf, graphml | Format of the output graph. |
| `--gff` | 0 | False | No | — | Generate GFF files with projected pangenome annotations for each input genome. |
| `--proksee` | 0 | False | No | — | Generate JSON map files for PROKSEE with projected pangenome annotations for each input genome. |
| `--table` | 0 | False | No | — | Generate a tsv file for each input genome with pangenome annotations. |
| `--compress` | 0 | False | No | — | Compress the files in .gz |
| `--add_sequences` | 0 | False | No | — | Include input genome DNA sequences in GFF and Proksee output. |
| `-c, --cpu` | int | 1 | No | — | Number of available cpus |
| `--tmpdir` | Path | /tmp | No | — | directory for storing temporary files |
| `--keep_tmp` | 0 | False | No | — | Keeping temporary files (useful for debugging). |
| `--add_metadata` | 0 | False | No | — | Include metadata information in the output files if any have been added to pangenome elements (see ppanggolin metadata command). |
| `--metadata_sources` | list | — | No | — | Which source of metadata should be written. By default all metadata sources are included. |
| `--metadata_sep` | str | `|` | No | — | The separator used to join multiple metadata values for elements with multiple metadata values from the same source. This character should not appear in metadata values. |
| `--verbose` | int | 1 | No | 0, 1, 2 | Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug) |
| `--log` | check_log | `stdout` | No | — | log output file |
| `-d, --disable_prog_bar` | 0 | False | No | — | disables the progress bars |
| `-f, --force` | 0 | False | No | — | Force writing in output directory and in pangenome output file. |
| `--config` | FileType('r') | — | No | — | Specify command arguments through a YAML configuration file. |

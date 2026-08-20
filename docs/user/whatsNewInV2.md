# What's new in PPanGGOLiN v2

PPanGGOLiN v2 is the result of the work carried out since mid-2021, released under version numbers `2.0.0` and above from January 2024 onwards. This page summarizes what changed with respect to PPanGGOLiN as originally published ([Gautreau et al. 2020](https://doi.org/10.1371/journal.pcbi.1007732)), taking release `1.1.136` (February 2021) as the reference point. It is intended for users coming from a v1 release who want to know what is new and what has moved.

The features listed below were released progressively, and some of them are already available in the later `1.2.x` releases. The major version number was raised to `2.0.0` at the point where the redesigned pangenome file broke compatibility with the HDF5 files produced by v1, marking the boundary between the two versions. The rationale behind these changes is described in the PPanGGOLiN v2 publication.

## New analyses

| Command       | Description                                                                                                                                                                                                                | Documentation                          |
|---------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------|
| `projection`  | Annotate external genomes using an already computed pangenome, without recomputing it. Genes are assigned to existing gene families and partitions, then RGPs, spots and modules are predicted for the projected genome.   | [Projection](projection.md)            |
| `context`     | Extract conserved genomic neighborhoods around genes of interest, using the neighborhood relationships already encoded in the pangenome graph.                                                                             | [Genomic context](genomicContext.md)   |
| `rgp_cluster` | Cluster RGPs across genomes based on their shared gene families, using a Gene Repertoire Relatedness (GRR) score.                                                                                                          | [RGP clustering](RGP/rgpClustering.md) |
| `metadata`    | Attach arbitrary metadata to any pangenome element (genes, contigs, genomes, gene families, edges, RGPs, spots and modules) from a tabulated file. Metadata is stored in the pangenome file and propagated to the outputs. | [Metadata](metadata.md)                |

## New utility commands

| Command          | Description                                                                        | Documentation                                                    |
|------------------|------------------------------------------------------------------------------------|------------------------------------------------------------------|
| `all`            | Run the complete pipeline, including RGPs, spots and modules, in a single command. | [Workflows](PangenomeAnalyses/pangenomeWorkflow.md)              |
| `metrics`        | Compute pangenome metrics such as genomic fluidity.                                | [Pangenome metrics](PangenomeAnalyses/pangenomeMetric.md)        |
| `write_metadata` | Export the metadata stored in a pangenome.                                         | [Metadata](metadata.md)                                          |
| `utils`          | Generate a default YAML configuration file with `--default_config`.                | [Configuration file](practicalInformation.md#configuration-file) |

## Reorganized commands

The `write` command has been split into three commands with clearer scopes:

| v1                                       | v2                                                                                                                           |
|------------------------------------------|------------------------------------------------------------------------------------------------------------------------------|
| `ppanggolin write --regions --spots ...` | `ppanggolin write_pangenome` — pangenome-level outputs (gene families, partitions, RGPs, spots, modules, graphs, statistics) |
| `ppanggolin write --projection`          | `ppanggolin write_genomes --table` — genome-level outputs, one file per genome                                               |
| —                                        | `ppanggolin write_metadata` — metadata outputs                                                                               |


## New output formats

- **Proksee** (`write_genomes --proksee`, and the `projection` outputs): JSON files that can be loaded in [Proksee](https://proksee.ca/) to explore circular genome maps annotated with pangenome information.
- **GFF3** (`write_genomes --gff`): genome annotations enriched with pangenome information (gene family, partition, RGP, spot, module).
- **graph-tool** (`write_pangenome --gt`): the pangenome graph in the compressed `gt` format, in addition to the GEXF and JSON formats already available in v1.
- **RGP gene families** (`write_pangenome --regions_families`): the gene family composition of each RGP.

## Usability

- **Configuration files.** Every command accepts a `--config` option pointing to a YAML file, so that parameters no longer have to be passed on the command line. This is particularly useful for the workflow commands, which run several steps in a single call, and for integration into computational platforms. A template can be generated with `ppanggolin utils --default_config`, and the format is described in the [Configuration file](practicalInformation.md#configuration-file) section.
- **Parameter tracking.** Parameters are resolved consistently (command line, then configuration file, then defaults), checked for consistency across analysis steps, and stored inside the pangenome file, so that an analysis can be inspected and reproduced afterward with `ppanggolin info --parameters`.
- **Defragmentation is now enabled by default.** The v1 opt-in `--defrag` option of the `cluster` command has been replaced by the opt-out `--no_defrag` option.
- **Documentation.** The documentation has been fully rewritten and is now hosted on [ReadTheDocs](https://ppanggolin.readthedocs.io), including a reference of all commands and their options and an API reference for programmatic use.
- **Distribution.** PPanGGOLiN v2 is available on [Bioconda](https://anaconda.org/bioconda/ppanggolin) and [PyPI](https://pypi.org/project/PPanGGOLiN/).

## Performance and file format

- **Pyrodigal replaces Prodigal.** Gene prediction now uses [Pyrodigal](https://github.com/althonos/pyrodigal), which removes the input/output overhead of calling an external binary during annotation.
- **Redesigned HDF5 structure.** The internal layout of the pangenome file was reworked, substantially reducing file size and giving faster and more convenient programmatic access.
- **Faster pangenome loading.** Reading a pangenome file has been optimized, reducing loading times from minutes to seconds on large datasets.
- **Lower memory usage when writing sequences.** Sequences are read directly from the HDF5 file instead of being loaded into memory beforehand.

## Compatibility with v1 pangenome files

The redesign of the HDF5 structure is what defines the boundary between v1 and v2: the major version number was raised to `2.0.0` precisely because pangenome files became incompatible at that point. Pangenome files produced with v1 therefore cannot be read by v2. PPanGGOLiN checks the version stored in the file and raises an explicit error if it was created by a v1 release, so such pangenomes have to be recomputed with v2.

## Citation

If you use features introduced in v2, please cite the PPanGGOLiN v2 publication in addition to the original PPanGGOLiN, panRGP and panModule papers listed on the [home page](../index.md).
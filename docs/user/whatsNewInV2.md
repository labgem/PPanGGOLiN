# What's new in PPanGGOLiN v2

This page summarises the changes between PPanGGOLiN v1 (last release `1.2.105`, February 2023) and v2 (first release `2.0.0`, January 2024). It is intended for users coming from v1 who want to know what is new and what has moved. The rationale behind these changes is described in the PPanGGOLiN v2 publication.

## New commands

| Command          | Description                                                                                                                                                                                                                | Documentation                                    |
|------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------|
| `projection`     | Annotate external genomes using an already computed pangenome, without recomputing it. Genes are assigned to existing gene families and partitions, then RGPs, spots and modules are predicted for the projected genome.   | [Projection](projection.md)                      |
| `rgp_cluster`    | Cluster RGPs across genomes based on their shared gene families, using a Gene Repertoire Relatedness (GRR) score.                                                                                                          | [RGP clustering](RGP/rgpClustering.md)           |
| `metadata`       | Attach arbitrary metadata to any pangenome element (genes, contigs, genomes, gene families, edges, RGPs, spots and modules) from a tabulated file. Metadata is stored in the pangenome file and propagated to the outputs. | [Metadata](metadata.md)                          |
| `write_metadata` | Export the metadata stored in a pangenome.                                                                                                                                                                                 | [Metadata](metadata.md)                          |
| `utils`          | Generate a default YAML configuration file with `--default_config`.                                                                                                                                                        | [Practical information](practicalInformation.md) |

## Reorganised commands

The v1 `write` command has been split into three commands with clearer scopes:

| v1                                                 | v2                                                                                                                           |
|----------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------|
| `ppanggolin write --regions --spots --modules ...` | `ppanggolin write_pangenome` — pangenome-level outputs (gene families, partitions, RGPs, spots, modules, graphs, statistics) |
| `ppanggolin write --projection`                    | `ppanggolin write_genomes --table` — genome-level outputs, one file per genome                                               |
| —                                                  | `ppanggolin write_metadata` — metadata outputs                                                                               |

Scripts written for v1 that call `ppanggolin write` need to be updated accordingly.

## New output formats

- **Proksee** (`write_genomes --proksee`, and the `projection` outputs): JSON files that can be loaded in [Proksee](https://proksee.ca/) to explore circular genome maps annotated with pangenome information.
- **GFF3** (`write_genomes --gff`): genome annotations enriched with pangenome information (gene family, partition, RGP, spot, module).
- **graph-tool** (`write_pangenome --gt`): the pangenome graph in the compressed `gt` format, in addition to the GEXF and JSON formats already available in v1.
- **RGP gene families** (`write_pangenome --regions_families`): the gene family composition of each RGP.

## Usability

- **Configuration files.** Every command accepts a `--config` option pointing to a YAML file, so that parameters no longer have to be passed on the command line. This is particularly useful for the workflow commands, which run several steps in a single call, and for integration into computational platforms. A template can be generated with `ppanggolin utils --default_config`.
- **Parameter tracking.** Parameters are resolved consistently (command line, then configuration file, then defaults), checked for consistency across analysis steps, and stored inside the pangenome file, so that an analysis can be inspected and reproduced afterwards with `ppanggolin info --parameters`.
- **Defragmentation is now enabled by default.** The v1 opt-in `--defrag` option of the `cluster` command has been replaced by the opt-out `--no_defrag` option.
- **Documentation.** The documentation has been fully rewritten and is now hosted on [ReadTheDocs](https://ppanggolin.readthedocs.io), including a reference of all commands and their options and an API reference for programmatic use.
- **Distribution.** PPanGGOLiN v2 is available on [Bioconda](https://anaconda.org/bioconda/ppanggolin) and [PyPI](https://pypi.org/project/PPanGGOLiN/).

## Performance and file format

- **Pyrodigal replaces Prodigal.** Gene prediction now uses [Pyrodigal](https://github.com/althonos/pyrodigal), which removes the input/output overhead of calling an external binary during annotation.
- **Redesigned HDF5 structure.** The internal layout of the pangenome file was reworked, substantially reducing file size and giving faster and more convenient programmatic access.
- **Faster pangenome loading.** Reading a pangenome file has been optimised, reducing loading times from minutes to seconds on large datasets.
- **Lower memory usage when writing sequences.** Sequences are read directly from the HDF5 file instead of being loaded into memory beforehand.

## Compatibility with v1 pangenome files

The HDF5 structure changed between v1 and v2. Pangenome files produced with v1 cannot be read by v2: PPanGGOLiN checks the version stored in the file and raises an explicit error if it was created by a v1 release. Such pangenomes have to be recomputed with v2.

## Citation

If you use features introduced in v2, please cite the PPanGGOLiN v2 publication in addition to the original PPanGGOLiN, panRGP and panModule papers listed on the [home page](../index.md).
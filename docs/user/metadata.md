table,# Metadata and Pangenome

## Associating Metadata to Pangenome Elements

The `metadata` command allows the addition of metadata linked to various pangenome elements. Metadata can be associated with genes, genomes, families, RGPs, spots, and modules using a simple TSV file. 

To add metadata to your pangenome, execute the command as shown below:

```bash
ppanggolin metadata -p PANGENOME --metadata METADATA.TSV --source SOURCE --assign ASSIGN
```

- The `--source` argument corresponds to the metadata's origin and will serve as the storage key in the pangenome.
- `--assign` allows you to specify the pangenome elements to which you want to add metadata from the following list: {families, genomes, genes, RGPs, spots, modules}.


The associated metadata can then be exported in various output files of PPanGGOLiN such as GFF, PROKSEE JSON Map and Table output for genomes (see [here](./writeGenomes.md#incorporating-metadata-into-tables-gff-and-proksee-files) for more details) and in the gexf graph file of the pangenome as well as in the graph resulting in the RGP clustering.


The metadata linked to pangenome elements can be exported to various output file formats within PPanGGOLiN, including GFF, PROKSEE JSON Map, and Table outputs of the `write_genomes` command (see [here](./writeGenomes.md#incorporating-metadata-into-tables-gff-and-proksee-files) for more details). Additionally, the metadata can also be included in the gexf graph file representing the pangenome and in the RGP clustering graph.

```{note}
Certain information (such as species, strain, and dbx_ref) is extracted from genome annotation files (GBFF, GFF) during the annotation step and added to the pangenome as metadata under the source 'annotation_files'. These metadata can be extracted using the `write_metadata` command.
```

### Metadata Format

PPanGGOLiN offers a highly flexible metadata file format. Only one column name is mandatory, and it aligns with the assignment argument chosen by the user (ie. families, RGPS...).

For instance, the TSV file used to assign metadata to gene families for functional annotation might resemble the following:

| families | Accession | Function | Description |
|----------|-----------|----------|-------------|
| GF_1     | Acc_1     | Fn_1     | Desc_1      |
| GF_2     | Acc_2     | Fn_2     | Desc_2      |
| GF_2     | Acc_3     | Fn_3     | Desc_3      |
| ...      | ...       | ...      | ...         |
| GF_n     | Acc_n     | Fn_n     | Desc_n      |

```{note} 
As you can see in the above table, one element (here GF_2) can be associated with with multiple metadata entries.
```

The file has to respect the following rules:

- it is tab-separated and has a header line;
- it contains at least two columns, one of them named after the `--assign` value;
- every column name is a valid identifier, that is, it only contains letters, digits and underscores, and does not start with a digit. Column names such as `E-value`, `hmm acc` or `target name` are rejected and have to be renamed;
- a cell containing a single `-` is interpreted as a missing value.

### Example: annotating gene families with Pfam

The following example shows the complete process of annotating the gene families of a pangenome with [Pfam](https://www.ebi.ac.uk/interpro/), and adding the result to the pangenome.

**1. Export the representative sequences of the gene families**

Each gene family of the pangenome is represented by one sequence, which can be written with the `fasta` command (see [how to write Gene families sequences](./writeFasta.md#gene-families) for more details):

```bash
ppanggolin fasta -p pangenome.h5 --output MY_PROT --prot_families all
```

This creates the file `MY_PROT/all_protein_families.faa`, containing one protein sequence per gene family. The identifier of each sequence is the identifier of the gene family it represents, which is precisely what has to be given back to PPanGGOLiN in the next steps.

**2. Annotate the sequences with Pfam**

The resulting FASTA file can then be given to any annotation tool. For Pfam, this is typically done with [HMMER](http://hmmer.org/) against the `Pfam-A.hmm` database, or with the `pfam_scan.pl` script distributed by Pfam:

```bash
hmmsearch --cut_ga --tblout pfam_hits.txt Pfam-A.hmm MY_PROT/all_protein_families.faa
```

**3. Format the result as a PPanGGOLiN metadata file**

The output of the annotation tool has to be turned into a TSV file respecting the rules given above. In practice this means keeping the columns of interest, renaming the column holding the gene family identifiers to `families`, and renaming any column whose name is not a valid identifier.

For instance, an annotation result looking like this:

| target_name | accession | query_name | E-value | score | description_of_target |
|-------------|-----------|------------|---------|-------|-----------------------|
| ABC_tran    | PF00005.30| GF_1       | 1.2e-45 | 152.3 | ABC transporter       |
| MFS_1       | PF07690.19| GF_2       | 3.4e-30 | 101.7 | Major facilitator superfamily |

becomes:

| families | pfam_name | pfam_accession | e_value | score  | description                   |
|----------|-----------|----------------|---------|--------|-------------------------------|
| GF_1     | ABC_tran  | PF00005.30     | 1.2e-45 | 152.3  | ABC transporter               |
| GF_2     | MFS_1     | PF07690.19     | 3.4e-30 | 101.7  | Major facilitator superfamily |

Only the column names change: the `query_name` column, which holds the gene family identifiers written at step 1, is renamed to `families`, and `E-value` is renamed to `e_value` so that it is a valid identifier. The remaining columns are free, both in number and in name.

```{note}
The order of the columns does not matter. PPanGGOLiN identifies the columns by their name, so the `families` column can be anywhere in the file. It has been moved to the first position in the example above only for readability.
```

**4. Add the annotations to the pangenome**

```bash
ppanggolin metadata -p pangenome.h5 --metadata pfam_families.tsv --source pfam --assign families
```

The annotations are now stored in the pangenome under the source `pfam`, and are propagated to the outputs described above. A gene family matched by several Pfam entries simply appears on several lines, as shown in the previous section.

### Command Specific Option Details

#### `--metadata`
PPanGGOLiN enables to give one TSV at a time to add metadata.

#### `--source`
The source serves as the key for accessing metadata within the pangenome. If the source name already exists in the pangenome, it can be overwritten using the `--force` option. This system facilitates the utilization of multiple metadata sources, accessible and usable within PPanGGOLiN. In the context of annotation, the source typically represents the name of the annotation database used during the annotation process. 

#### `--assign`
PPanGGOLiN enables the addition of metadata to various pangenome elements, including families, genomes, genes, RGPs, spots, and modules. However, the user can provide only one metadata file at a time, thereby specifying a single source and one type of pangenome element.

#### `--omit`
This option allows you to bypass errors resulting from an unfound ID in the pangenome. It can be beneficial when utilizing a general TSV with elements not present in the pangenome. This argument should be used carefully.  

## Checking Metadata Associated with the Pangenome

You can inspect which pangenome elements have associated metadata and their respective sources using the `ppanggolin info` command. For more details on this command, refer to the [info command documentation](./PangenomeAnalyses/pangenomeInfo.md).

To check the metadata, execute the following command:

```bash
ppanggolin info -p PANGENOME --metadata
```

## Exporting Metadata to TSV Files

The `write_metadata` command allows you to export metadata to TSV files. It creates one file per source of pangenome element metadata.

For example, if families have metadata from sources `hmmer` and `dbcan`, and genomes have metadata from the source `gtdb`, the command will create three files:
- `families_metadata_from_hmmer.tsv`
- `families_metadata_from_dbcan.tsv`
- `genomes_metadata_from_gtdb.tsv`

To export metadata from your pangenome, execute the following command:

```bash
ppanggolin write_metadata -p PANGENOME --output metadata_tsv_files
```

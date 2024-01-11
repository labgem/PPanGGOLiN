# Adding Metadata to Pangenome Elements

The `metadata` command allows the addition of metadata linked to various pangenome elements. Metadata can be associated with genes, genomes, families, RGPs, spots, and modules using a simple TSV file. 

To add metadata to your pangenome, execute the command as shown below:

```bash
ppanggolin metadata -p PANGENOME --metadata METADATA.TSV --source SOURCE --assign ASSIGN
```

- The `--source` argument corresponds to the metadata's origin and will serve as the storage key in the pangenome.
- `--assign` allows you to specify the pangenome elements to which you want to add metadata from the following list: {families, genomes, genes, RGPs, spots, modules}.


The associated metadata can then be exported in various output files of PPanGGOLiN such as GFF, PROKSEE JSON Map and Table output for genomes (see [here](./writeGenomes.md#incorporating-metadata-into-tables-gff-and-proksee-files) for more details) and in the gexf graph file of the pangenome as well as in the graph resulting in the RGP clustering.


The metadata linked to pangenome elements can be exported to various output file formats within PPanGGOLiN, including GFF, PROKSEE JSON Map, and Table outputs of the `write_genomes` command (see [here](./writeGenomes.md#incorporating-metadata-into-tables-gff-and-proksee-files) for more details). Additionally, the metadata can also be included in the gexf graph file representing the pangenome and in the RGP clustering graph.

## Metadata Format

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

### Command Specific Option Details

#### `--metadata`
PPanGGOLiN enables to give one TSV at a time to add metadata.

#### `--source`
The source serves as the key for accessing metadata within the pangenome. If the source name already exists in the pangenome, it can be overwritten using the `--force` option. This system facilitates the utilization of multiple metadata sources, accessible and usable within PPanGGOLiN. In the context of annotation, the source typically represents the name of the annotation database used during the annotation process. 

#### `--assign`
PPanGGOLiN enables the addition of metadata to various pangenome elements, including families, genomes, genes, RGPs, spots, and modules. However, the user can provide only one metadata file at a time, thereby specifying a single source and one type of pangenome element.

#### `--omit`
This option allows you to bypass errors resulting from an unfound ID in the pangenome. It can be beneficial when utilizing a general TSV with elements not present in the pangenome. This argument should be used carefully.  

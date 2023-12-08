# Add metadata to pangenome elements

It is possible to add metadata link to pangenome elements using PPanGGOLiN. 
Metadata can be associated with: genes, genomes, families, RGPs, spots and modules from a simple TSV file. 
To add metadata in your pangenome you can launch the command is as follows:

`ppanggolin metadata -p PANGENOME --metadata METADATA.TSV --source SOURCE --assign ASSIGN`

- `--source` arguments corresponds to the origin of the metadata and will be used as the storage key in the pangenome.
- `--assign` Choose to which pangenome elements who want to add metadata in the following list {families,genomes,genes,RGPs,spots,modules}

## Metadata format

PPanGGOLiN allows to use a highly flexible metadata file. Only one column name is mandatory, and it is identical to the 
assignment argument chosen by the user.

For example the TSV file to assign metadata to gene families to functional annotation could be as follows:

| families | Accesion | Function | Description |
|----------|----------|----------|-------------|
| GF_1     | Acc_1    | Fn_1     | Desc_1      |
| GF_2     | Acc_2    | Fn_2     | Desc_2      |
| GF_2     | Acc_3    | Fn_3     | Desc_3      |
| ...      | ...      | ...      | ...         |
| GF_n     | Acc_n    | Fn_n     | Desc_n      |

*Note: As you can see in the above table, one element (here GF_2) can be associated with more than one metadata.*

### Command specifiq option details

#### `--metadata`
PPanGGOLiN enables to give one TSV at a time to add metadata.

#### `--source` 
The source is the key use to access to metadata in pangenome. 
So if the name of the source already exist in the pangenome it can be overwritten only with `--force` option.
This system allow to have multiple metadata source that can be read and use in PPanGGOLiN.

#### `--assign` 
PPanGGOLiN allows to add metadata to all pangenome elements: families,genomes,genes,RGPs,spots,modules.
But the user can only give one metadata file at a time as he can provide only source and so one type of pangenome element.

#### `--omit`
You can use this option to skip the error provide by an unfind ID in the pangenome. 
This could be useful if you are using a general TSV with element not in the pangenome, but must be used with carefully.  


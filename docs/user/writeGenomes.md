# Write genomes

The `write_genomes` command creates 'flat' files representing genomes with their pangenome annotations.

To generate output for specific genomes, use the `--organisms` argument. This argument accepts a list of organism names, either directly entered in the command line (comma-separated) or referenced from a file where each line contains a single organism name.


### Genes table with pangenome annotations

The `--table` option generates a TSV file for each genome, which contains pangenome annotations for the genes.  These files are stored in a directory named `tables`.


The following table outlines the columns present in the generated files:

| Column               | Description                                                                |
|----------------------|----------------------------------------------------------------------------|
| gene                 | Unique identifier for the gene                                             |
| contig               | Contig on which the gene is located                                        |
| start                | Start position of the gene                                                 |
| stop                 | Stop position of the gene                                                  |
| strand               | Gene location strand                                      |
| family               | ID of the gene's associated family in the pangenome             |
| nb_copy_in_org       | Number of copies of a family present in the organism; 1 indicates no close paralogs |
| partition            | Gene family partition in the pangenome                  |
| persistent_neighbors | Number of neighbors classified as 'persistent' in the pangenome graph        |
| shell_neighbors      | Number of neighbors classified as 'shell' in the pangenome graph             |
| cloud_neighbors      | Number of neighbors classified as 'cloud' in the pangenome graph             |
| RGP                  | ID of the Region of Genomic Plasticity (RGP) if the gene is within an RGP <br> (present only if RGPs have been predicted) |
| spot                 | Spot ID in which the RGP is inserted (present only if RGPs and spot have been predicted)   |
| module               | Module ID of the gene family (present if modules have been predicted)             |


```{note}
Columns such as RGP, spot, and module are included only when these elements have been predicted in the pangenome.
```

Those files can be generated as such : 

```bash
ppanggolin write_genomes -p pangenome.h5 --table -o write_genomes_output
```


### GFF file


The `--gff` argument generates GFF files, each containing pangenome annotations for individual genomes within the pangenome. The GFF file format is a widely recognized standard in bioinformatics and can easily integrate into downstream analysis tools.

To generate GFF files from a pangenome HDF5 file, you can use the following command:

```bash
ppanggolin write_genomes -p pangenome.h5 --gff -o output
```

This command will create a gff directory within the output directory, with one GFF file per genome. 

Pangenome annotations within the GFF are recorded in the attribute column of the file.

For CDS features, pangenome annotations are recorded in the attribute column of the file:

CDS features have the following attributes:

- **family:** ID of the gene family to which the gene belongs.
- **partition:** Partition of the gene family, categorized as persistent, shell, or cloud.
- **module:** If the gene family belongs to a module, the module ID is specified with the key 'module.'
- **rgp:** If the gene is part of a Region of Genomic Plasticity (RGP), the RGP name is specified with the key 'rgp.'

For Regions of Genomic Plasticity (RGPs), RGPs are specified under the feature type 'region.'

RGPs have the following attributes:

- **'spot'** designates the spot ID where the RGP is inserted. When the RGP has no spot, the term 'No_spot' is used.
- **'Note'** specifies that this feature is an RGP.


Here is an example showcasing the initial lines of the GFF file for the _Acinetobacter baumannii_ AYE genome:

```gff
##gff-version 3
##sequence-region NC_010401.1 1 5644
##sequence-region NC_010402.1 1 9661
##sequence-region NC_010403.1 1 2726
##sequence-region NC_010404.1 1 94413
##sequence-region NC_010410.1 1 3936291
NC_010401.1	.	region	1	5644	.	+	.	ID=NC_010401.1;Is_circular=true
NC_010401.1	ppanggolin	region	629	5591	.	.	.	Name=NC_010401.1_RGP_0;spot=No_spot;Note=Region of Genomic Plasticity (RGP)
NC_010401.1	external	gene	629	1579	.	+	.	ID=gene-ABAYE_RS00005
NC_010401.1	external	CDS	629	1579	.	+	0	ID=ABAYE_RS00005;Parent=gene-ABAYE_RS00005;product=replication initiation protein;family=ABAYE_RS00005;partition=cloud;rgp=NC_010401.1_RGP_0
NC_010401.1	external	gene	1576	1863	.	+	.	ID=gene-ABAYE_RS00010
NC_010401.1	external	CDS	1576	1863	.	+	0	ID=ABAYE_RS00010;Parent=gene-ABAYE_RS00010;product=hypothetical protein;family=ABAYE_RS00010;partition=cloud;rgp=NC_010401.1_RGP_0
NC_010401.1	external	gene	2054	2572	.	-	.	ID=gene-ABAYE_RS00015
NC_010401.1	external	CDS	2054	2572	.	-	0	ID=ABAYE_RS00015;Parent=gene-ABAYE_RS00015;product=tetratricopeptide repeat protein;family=HTZ92_RS18670;partition=shell;rgp=NC_010401.1_RGP_0
```

### JSON Map for Proksee visualisation


The `--proksee` argument generates JSON map files containing pangenome annotations, which can be visualized using Proksee at [https://proksee.ca/](https://proksee.ca/).

To generate JSON map files, you can use the following command:

```bash
ppanggolin write_genomes -p pangenome.h5 --proksee -o output
```

This command creates a proksee directory within the output directory and generates one JSON file per genome. 


To load a JSON map file on Proksee, follow these steps:
1. Navigate to the "Map JSON" tab.
2. Upload your file using the browse button.
3. Click the "Create Map" button to generate the visualization.

A genome visualized by Proksee with PPanGGOLiN annotation appears as depicted below:


```{image} ../_static/proksee_exemple_A_baumannii_AYE.png
:align: center
```

*Image: Genome visualized by Proksee with PPanGGOLiN annotation.*


The visualization is composed of three tracks:
- **Genes:** Color-coded by their gene family partition.
- **RGP:** The annotation of the object specifies the spot associated with the RGPs.
- **Module:** Displaying modules within the genome. The completion of the module is specified in the annotation of the object.





### Adding Fasta Sequences into GFF and proksee JSON map Files

PPanGGOLiN allows the incorporation of fasta sequences into GFF files and proksee JSON map files. This integration with Proksee provides access to various tools that rely on DNA sequences, including the construction of GC% and GC skew profiles, and conducting blast searches for example.


Since PPanGGOLiN does not retain genomic sequences, it is necessary to provide the original genomic files used to construct the pangenome through either the `--anno` or `--fasta` argument. These arguments mirror those used in workflow commands (`workflow`, `all`, `panrgp`, `panmodule`) and the `annotate` command.

- `--anno`: This option requires a tab-separated file containing organism names and the corresponding GFF/GBFF file paths of their annotations. If `--anno` is utilized, GFF files should include fasta sequences.

- `--fasta`: Use this option with a tab-separated file that lists organism names alongside the filepaths of their genomic sequences in fasta format.


### Incorporating Metadata into Tables, GFF, and Proksee Files

You can inject metadata, previously added with the `metadata` command, into genome outputs using the `--add_metadata` parameter. When users add metadata, they specify the source of this metadata. These metadata sources can be selectively included using the `--metadata_sources` parameter. By default, all sources are added when the `--add_metadata` flag is specified.

#### Metadata in GFF Files

Metadata is integrated into the attributes column of the GFF file. The patterns for adding metadata are as follows:

- In CDS lines, metadata associated with genes follow this pattern: `gene_<source>_<key>=<value>`. Gene family metadata follows a similar pattern: `gene_<source>_<key>=<value>`.
- In the contig lines of type `region` describing the contig, genome metadata is added with the pattern: `genome_<source>_<key>=<value>`, and contig metadata is added with: `contig_<source>_<key>=<value>`.
- In RGP lines, metadata is added using the pattern: `rpg_<source>_<key>=<value>`.

For example, if we associate the metadata of the gene family DYB08_RS16060 with the source `pfam`:

```tsv
families	accession	type	description
DYB08_RS16060	PF18894	domain	This entry represents a probable metallopeptidase domain found in a variety of phage and bacterial proteomes.
```

This metadata file can be added to the pangenome with the metadata command:

```bash
ppanggolin metadata -p pangenome.h5 --source pfam --metadata family_pfam_annotation.tsv --assign families
```

When writing GFF output with the `--add_metadata` flag:

```bash
ppanggolin write_genomes -p pangenome.h5 --proksee -o proksee_out --gff --add_metadata
```

A gene belonging to this family would have the following attribute in its GFF line: `family_pfam_accession=PF18894;family_pfam_description=This entry represents a probable metallopeptidase domain found in a variety of phage and bacterial proteomes.;family_pfam_type=domain`.

```gff
NC_010404.1	external	CDS	77317	77958	.	-	0	ID=ABAYE_RS00475;Parent=gene-ABAYE_RS00475;product=putative metallopeptidase;family=DYB08_RS16060;partition=persistent;rgp=NC_010404.1_RGP_0;family_pfam_accession=PF18894;family_pfam_description=This entry represents a probable metallopeptidase domain found in a variety of phage and bacterial proteomes.;family_pfam_type=domain
```

#### Metadata in Proksee Visualization

Metadata can be seamlessly incorporated into Proksee JSON MAP files. These metadata details become accessible by simply hovering the mouse over the features.

For instance, with the metadata previously added to the DYB08_RS16060 gene family, the Proksee visualization would look like the following figure:

```{image} ../_static/proksee_metadata_example.png
:align: center
```


#### Metadata in Table output

Metadata is seamlessly incorporated into table output with the addition of extra columns. These columns follow the GFF attribute naming: 

- gene metadata: `gene_<source>_<key>`
- family metadata: `gene_<source>_<key>`

<!-- exemple ? -->


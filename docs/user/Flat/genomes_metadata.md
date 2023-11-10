<!-- ### Incorporating Metadata into Tables, GFF, and Proksee Files -->

You can inject metadata, previously added with the `metadata` command, into genome outputs using the `--add_metadata` parameter. When users add metadata, they specify the source of this metadata. These metadata sources can be selectively included using the `--metadata_sources` parameter. By default, all sources are added when the `--add_metadata` flag is specified.

#### Metadata in GFF Files

Metadata is integrated into the attributes column of the GFF file. The patterns for adding metadata are as follows:

- In CDS lines, metadata associated with genes follow this pattern: `gene_<source>_<column>=<value>`. Gene family metadata follows a similar pattern: `gene_<source>_<column>=<value>`.
- In the contig lines of type `region` describing the contig, genome metadata is added with the pattern: `genome_<source>_<column>=<value>`, and contig metadata is added with: `contig_<source>_<column>=<value>`.
- In RGP lines, metadata is added using the pattern: `rpg_<source>_<column>=<value>`.

For example, if we associate metadata is associated with the gene family DYB08_RS16060 with the source `pfam`:

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

### Metadata in Proksee Visualization

Metadata can be seamlessly incorporated into Proksee JSON MAP files, enriching the visualization experience. These metadata details become accessible by simply hovering the mouse over the features.

For instance, with the metadata previously added to the DYB08_RS16060 gene family, the Proksee visualization would resemble the example below:

```{image} ../_static/proksee_metadata_example.png
:align: center
```

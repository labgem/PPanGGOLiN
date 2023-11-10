
The `--gff` argument generates GFF files, each containing pangenome annotations for individual genomes within the pangenome. The GFF file format is a widely recognized standard in bioinformatics and can seamlessly integrate into downstream analysis tools.

To generate GFF files from a pangenome HDF5 file, you can use the following command:

```bash
ppanggolin write_genomes -p pangenome.h5 --gff -o output
```

This command will create a gff directory within the output directory, with one GFF file per genome. 

Pangenome annotations within the GFF are recorded in the attribute column of the file.

For CDS features, pangenome annotations are recorded in the attribute column of the file:

CDS features have the following attributes:

- **family:** ID of the gene family to which the gene belongs.
- **partition:** The partition of the gene family, categorized as persistent, shell, or cloud.
- **module:** If the gene family belongs to a module, the module ID is specified with the key 'module.'
- **rgp:** If the gene is part of a Region of Genomic Plasticity (RGP), the RGP name is specified with the key 'rgp.'

For Regions of Genomic Plasticity (RGPs), RGPs are specified under the feature type 'region.'

RGPs have the following attributes:

- The attribute 'spot' designates the spot ID where the RGP is inserted. When the RGP has no spot, the term 'No_spot' is used.
- The 'Note' attribute specifies that this feature is an RGP.


Here is an example showcasing the initial lines of the GFF file for the Acinetobacter baumannii AYE genome:

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
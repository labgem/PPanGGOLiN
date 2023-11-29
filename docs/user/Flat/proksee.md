The `--proksee` argument generates JSON map files containing pangenome annotations, which can be visualized using Proksee at [https://proksee.ca/](https://proksee.ca/).

To generate JSON map files, you can use the following command:

```bash
ppanggolin write_genomes -p pangenome.h5 --proksee -o output
```

This command will create a proksee directory within the output directory, with one JSON file per genome. 


To load a JSON map file on Proksee, follow these steps:
1. Navigate to the "Map JSON" tab.
2. Upload your file using the browse button.
3. Click the "Create Map" button to generate the visualization.

A genome visualized by Proksee with PPanGGOLiN annotation appears as depicted below:


```{image} ../_static/proksee_exemple_A_baumannii_AYE.png
:align: center
```

*Image: Genome visualized by Proksee with PPanGGOLiN annotation.*


The visualization consists of three tracks:
- **Genes:** Color-coded by their gene family partition.
- **RGP (Region of Genomic Plasticity):** Spot associated to the RGPs are specified in the annotation of the object.
- **Module:** Displaying modules within the genome. The completion of the module is specified in the annotation of the object.


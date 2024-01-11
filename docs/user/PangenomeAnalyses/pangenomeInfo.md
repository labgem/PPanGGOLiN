
<!-- ### Explore the Contents of Your Pangenome File -->

The `info` command in PPanGGOLiN enables users to acquire comprehensive insights into the contents and construction process of a pangenome file.

Different types of information can be displayed using various parameters, such as `--status`, `--parameters`, `--content`, and `--metadata`. When no flag is specified, all available outputs are displayed.

```bash
ppanggolin info -p pangenome.h5
```

The `info` command generates information in YAML format, displaying the results directly in the standard output without writing to any files.

#### Overview of `info --content` Output

The `--content` option with the `info` command exhibits general statistical data about your pangenome:

1. **General Metrics**:
   - Presents a comprehensive count of genes, genomes, gene families, and edges within the pangenome graph.

2. **Partitioned Metrics**:
   - Provides detailed information for each partition, including gene counts and presence thresholds for persistent, shell, and cloud families among genomes.

3. **Number of Partitions**:
   - Indicates the total count of partitions in the pangenome.

4. **Genome Fluidity**:
   - If computed through the 'metrics' command, displays genome fluidity values across all partitions and within each partition.

5. **Regions of Genomic Plasticity (RGPs)**:
   - Exhibits counts of Regions of Genomic Plasticity (RGPs) and spots of insertion if predicted using commands such as 'all', 'panrgp', 'rgp', or 'spot'.

6. **Modules**:
   - Shows counts of modules and associated gene families if predicted using commands like 'module', 'panmodule', or 'all'. Additionally, provides partition composition percentages for these modules.

#### Overview of `info --parameters` Output

This option displays the PPanGGOLiN parameters used at each analysis step. The output can be utilized as a configuration file for other PPanGGOLiN commands to replicate the same parameters. Refer [here](../practicalInformation.md#configuration-file) for more details on the configuration file .

#### Overview of `info --status` Output

Using this option, users can check what analysis have been conducted to obtain the pangenome file.

#### Overview of `info --metadata` Output

When metadata has been added to the pangenome elements, this option showcases which elements possess metadata and their respective sources. Find more details on metadata [here](../metadata.md).

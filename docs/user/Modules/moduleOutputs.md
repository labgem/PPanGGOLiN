## Module outputs


### Descriptive Tables for Predicted Modules

To describe predicted modules, various files can be generated, each describing different characteristics of these modules.

To generate these tables, use the `write_pangenome` command with the `--module` :

```bash
ppanggolin write_pangenome -p pangenome.h5 --modules -o my_output_dir
```

This command generates three tables: `functional_modules.tsv`, `modules_in_genomes.tsv`, and `modules_summary.tsv` described below:


#### 1. Gene Family to Module Mapping Table

The `functional_modules.tsv` file lists modules with their corresponding gene families. Each line establishes a mapping between a gene family and its respective module. 

It follows the following format:

|Column|Description|
|------|------------|
|module_id| Identifier for the module|
|family_id| Identifier for the family|


#### 2. Genome-wise Module Composition

The `modules_in_genomes.tsv` file provides a comprehensive overview of the modules present in each genome, detailing their completeness levels. Due to potential variability in module predictions, some modules might exhibit partial completeness in specific genomes where they are detected.

The structure of the `modules_in_genomes.tsv` file is outlined as follows:

| Column       | Description                                   |
|--------------|-----------------------------------------------|
| module_id    | Identifier for the module                      |
| genome       | Genome in which the indicated module is found  |
| completion   | Indicates the level of completeness (0.0 to 1.0) of the module in the <br> specified genome based on gene family representation |


#### 3. modules summary

The `modules_summary.tsv` file lists  characteristics for each detected module, with one line for each module.
The format is as follows:

|column|description|
|------|------------|
|module_id| The module identifier|
|nb_families| The number of families which are included in the module The families <br> themselves are listed in the 'functional_modules.tsv' file.|
|nb_genomes|The number of genomes in which the module is found. Those genomes are <br> listed in the 'modules_in_genomes.tsv' file.|
|partition| The average partition of the families in the module.|
|mean_number_of_occurrence| the mean number of time a module is present in each genome. <br> The expected value is around one, but it can be more if it is a module often repeated in the genomes (like a phage).|


### Mapping Modules with Spots and Regions of Genomic Plasticity (RGPs)

Predicted modules can be associated with Spots of insertion and Regions of Genomic Plasticity (RGPs) using the `write_pangenome` command with the `--spot_modules` flag as follows:

```bash
ppanggolin write_pangenome -p pangenome.h5 --spot_modules -o my_output_dir
```

This command generates two tables: `modules_spots.tsv` and `modules_RGP_lists.tsv`, described below.

```{note}
These outputs are available only if modules, spots, and RGPs have been computed in your pangenome (see the command [`all`](../QuickUsage/quickWorkflow.md#ppanggolin-complete-workflow-analyses) or the commands [`spot`](../RGP/rgpPrediction.md#spot-prediction), [`rgp`](../RGP/rgpPrediction.md#rgp-detection), and [`module`](./modulePrediction.md#conserved-module-prediction) for that).
```

Moreover, this information can be visualized through figures using the command `ppanggolin draw --spots` (refer to [Spot plots](../RGP/rgpOutputs.md#draw-spots), which can display modules).

#### 1. Associating Modules and Spots

The `modules_spots.tsv` file indicates which modules are present in each spot.

Its format is as follows:

| Column     | Description        |
|------------|--------------------|
| module_id  | Module identifier  |
| spot_id    | Spot identifier    |

#### 2. Associating Modules and RGPs

The `modules_RGP_lists.tsv` file lists RGPs that contain the same modules. These RGPs may have different gene families, but they will not include any other modules apart from those indicated. The format of `modules_RGP_lists.tsv` is as follows:

| Column             | Description                                                                                       |
|--------------------|---------------------------------------------------------------------------------------------------|
| representative_RGP | An RGP considered representative for the group, serving as a randomly chosen 'group of RGP IDs'   |
| nb_spots           | The number of spots where the RGPs containing the listed modules are observed                     |
| mod_list           | A list of the modules present in the indicated RGPs                                                 |
| RGP_list           | A list of RGPs that specifically include the previously listed modules                             |



### Module Information

To gather additional insights into the modules, including information about families and their distribution across different partitions, you can use the following command:

```bash
ppanggolin metrics -p pangenome.h5 --info_modules
```

The command output provides the following details:

```yaml
- Modules: 3
	- Families in Modules: 15
	- Percent of Families: 
		- persistent: 0.0
		- shell 53.33
		- cloud 46.67
	- Number of Families per Modules:
		- min: 3
		- max: 8
		- sd: 2.65
		- mean: 5
```


#### Functional modules
This .tsv file lists the modules and the gene families that belong to them. It lists one family per line, and there are multiple line for each module.
It is written along with other files with the following command:
`ppanggolin write -p pangenome.h5 --modules`

It follows the following format:
|column|description|
|------|------------|
|module_id| The module identifier|
|family_id| the family identifier|

#### Modules in organisms
This .tsv file lists for each organism the modules that are present and how complete they are. Since there are some variability that are allowed in the module predictions, occasionnally some modules can be incomplete in some of the organisms where they are found.
This file is written along with other files with the following command:
`ppanggolin write -p pangenome.h5 --modules`

And it follows the following format:
|column|description|
|------|------------|
|module_id| The module identifier|
|organism| the organism which has the indicated module|
|completion| a value between 0.0 and 1.0 which indicates how complete (in terms of gene family) the module is in the given organism|

#### modules summary
This .tsv file lists a few characteristics for each detected module. There is one line for each module.
The file is written along with other files with the following command:
`ppanggolin write -p pangenome.h5 --modules`

And it follows the following format:
|column|description|
|------|------------|
|module_id| The module identifier|
|nb_families| The number of families which are included in the module The families themselves are listed in the 'functional_modules.tsv' file.|
|nb_organisms|The number of organisms in which the module is found. Those organisms are listed in the 'modules_in_organisms.tsv' file.|
|partition| The average partition of the families in the module.|
|mean_number_of_occurrence| the mean number of time a module is present in each organism. The expected value is around one, but it can be more if it is a module often repeated in the genomes (like a phage).|

### spot modules
This command is available only if both modules and spots have been computed for your pangenome (see the command `all`, or the commands `spot` and `module` for that).
It indicates which modules are present in which spot and in which RGP.
The files are written with the following command:
```ppanggolin write -p pangenome.h5 --spot_modules```
The format of the 'modules_spots.tsv' file is the following:

|column|description|
|------|------------|
|module_id| The module identifier|
|spot_id| the spot identifier|

The file 'modules_RGP_lists.tsv' lists RGPs that have the same modules. Those RGPs can have different gene families, however they will not have any other module than those that are indicated. The format of the 'modules_RGP_lists.tsv' is the following:

|column|description|
|------|------------|
|representative_RGP| an RGP deemed representative for the group, and serving as a 'group of rgp id'(randomly picked)|
|nb_spots| The number of spots in which we see the RGPs which have the modules listed afterwards|
|mod_list| a list of the modules that are in the indicated RGPs|
|RGP_list| a list of RGP that include exactly the modules listed previously|

This information can also be visualized through figures that can be drawn with `ppanggolin draw --spots` (see [Spot plots](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#spot-plots), and which can display modules.

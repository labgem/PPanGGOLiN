# How to prepare your working environment

In order to work with PPanGGOLiN and to follow the tutorial we recommend to follow the next steps.

## Create a conda environment

The first step consist of creating a conda environment to install PPanGGOLiN and its dependencies.
You can look [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to know how install conda on your system. 

You can create tour environment as such:
```shell
conda create --name ppanggo 
```

More information on how to create a conda environment [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment)


## Install PPanGGOLiN

To install PPanGGOLiN in your conda environment you need to add 3 channels as follows:

```
conda config --add channels defaults ;
conda config --add channels bioconda ;
conda config --add channels conda-forge ;
```

Then you can just run:

```
conda install -c bioconda ppanggolin
```

This command will automatically install PPanGGOLiN dependencies. 
If you want more information or other methods to install PPanGGOLiN look at the documentation [here](../user/Installation.md#Installation)

## Install tutorial dependencies

As part of the tutorial we will also install some other software and packages.

### Download genomes

To download our genomes, we are going to use [genome_updater](https://github.com/pirovc/genome_updater).
Other solution exist such as [ncbi genome downloading scripts](https://github.com/kblin/ncbi-genome-download). Feel free to use the best and easiest way for you.

To install genome updater you can use the following command: `conda install -c bioconda genome_updater`

To install ncbi genome downloader you can use the following command: `conda install -c bioconda ncbi-genome-download`

### Graph visualisation

To visualise the pangenome graph you can use [Gephi](https://gephi.org/). To install it you must download the archive [here](https://gephi.org/users/download/). Then you can follow the next command:

```shell
tar -xvzf path/to/gephi/archive/gephi-X.XX.X-linux-x64.tar.gz
cd gephi-X.XX.X
chmod 755 bin/gephi
./bin/gephi
```

```{tip}
Gephi is also available on flathub and snapcraft.
```
[//]: # (Other depencies can be added here)
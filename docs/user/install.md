# Installation

## Latest version 

### Install with conda (recommended)

The recommended way of installing PPanGGOLiN is to do it using conda. 
To prevent conflicting dependencies, we also recommend to install PPanGGOLiN in its own environment as such

```
# Install into a new conda environment
conda create -n ppanggo -c conda-forge -c bioconda ppanggolin

# Check PPanGGOLiN install
conda activate ppanggo
ppanggolin --version
```

```{tip}
You can also use mamba, which is much quicker and sometimes help solving conflicting dependencies
```

Sometimes, installation problems come from having an unsupported python version installed by default in your environment.
Forcing the python version solves the problem:

``` 
conda create -n ppanggo -c conda-forge -c bioconda ppanggolin python=3.10
```

```{note}
Supported python version are 3.8, 3.9 and 3.10
```

### Install from source code (GitHub)

If you want to install from the source code, you should respect some step before.

```{warning}
This is a manual installation, we can help, but not support all troubleshooting. Prefer to use the conda installation whenever possible.
```

First, you must install the python dependencies.
For that, create a *dependencies.txt* file and copy/paste the next content.

```text
tqdm>=4.66.0
tables>=3.8.0
pyrodigal>=3.0.1
networkx>=3.0
scipy>=1.10.0
plotly>=5.16.0
gmpy2>=2.1.0
pandas>=2.0.0
colorlover>=0.3
numpy>=1.24.0
bokeh>=3.1.0
```

Next, you can use **pip** to install dependencies

```bash
python3 -m pip install -r dependencies.txt
```
```{warning}
Be sure to use a python version greater than 3.8 !
```
Then you must install the following software:

- [Aragorn>=1.2.41](http://www.ansikte.se/ARAGORN/Downloads/)
- [Infernal>=1.1.4](http://eddylab.org/infernal/)
- [MMSeqs2>=13.45111](https://github.com/soedinglab/MMseqs2/wiki#installation)

To finish, you can install ppanggolin by cloning the GitHub repository and using **pip**

```bash
git clone https://github.com/labgem/PPanGGOLiN.git
cd PPanGGOLiN
pip install .
```

## Development version

If you want to use the development version, you can use the 'dev' branch on GitHub.
While it is not guaranteed to work, it should most of the time.


You need to clone the repository on your computer, as followed:

```bash
git clone https://github.com/labgem/PPanGGOLiN.git
cd PPanGGOLiN
git checkout -b dev
git branch --set-upstream-to=origin/dev dev
git pull 
```

Then you will need to have all the dependencies.
You can install them as described [above](#install-from-source-code-github) with **pip** or by using **conda** and the [*requirements.txt file*](https://github.com/labgem/PPanGGOLiN/blob/dev/requirements.txt) as such:

```
conda create -n ppanggo -c conda-forge -c bioconda --file requirements.txt 
```

Then you will need to install ppanggolin. 
You can use pip while in your conda environment, which is the preferred way of doing. 
For this, place yourself at the root of the repository and run : 

```
pip install .
```

It should install PPanGGOLiN in your conda environment.
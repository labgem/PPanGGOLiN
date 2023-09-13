# Installation

## Latest version 

### Install with conda (recommended)

The recommended way of installing PPanGGOLiN is to do it using conda.To use it you need to add the conda channels that store the dependencies as such :

```
conda config --add channels defaults ;
conda config --add channels bioconda ;
conda config --add channels conda-forge ;
```

Then you can just run : 

`conda install -c bioconda ppanggolin`

You can also use mamba, which is much quicker and sometimes help solving conflicting dependencies:

```
conda install mamba
mamba install -c bioconda ppanggolin
```

If you have troubles or if conda tells you something about conflicting dependencies, I recommend you install PPanGGOLiN on a separate conda environment as PPanGGOLiN has quite a few dependencies and their versions can be conflicting with other bioinformatics software.

Sometimes, installation problems come from having an unsupported python version installed by default in your environment. Forcing the python version solves the problem:
``` 
conda install -c bioconda ppanggolin python=3.8
```

### Install from source code (GitHub)

If you want to install from the source code, you should respect some step before.

```{note}
This is a manual installation, we can help, but not support all troubleshooting. Prefer to use the conda installation whenever possible.  
```

First, you must install the python dependencies.
For that, create a *dependencies.txt* file and copy/paste the next content.

```text
tqdm>=4.66.0
tables>=3.8.0
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

- [Prodigal>=2.6.3](https://github.com/hyattpd/Prodigal/wiki/installation) 
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

If you want to use the development version you can use the 'dev' branch on GitHub. While it is not guaranteed to work, it should most of the time.

You need to clone the repository on your computer, as followed:

```bash
git clone https://github.com/labgem/PPanGGOLiN.git
cd PPanGGOLiN
git fetch
git checkout -b dev
git pull
```
then you will need to have all the dependencies.
You can install them as described above with **pip** or by using conda and the *requirements.txt file* as such:

` conda install --file requirements.txt `

Then you will need to install ppanggolin. You can use pip while in your conda environment which is the preferred way of doing. For this, place yourself at the root of the repository and run : 

`pip install .`

It should install PPanGGOLiN in your conda environment.
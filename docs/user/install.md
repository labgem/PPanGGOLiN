# Installation

## Installing PPanGGOLiN with Conda (recomended)

The best way to install PPanGGOLiN by utilizing conda, available through the bioconda channel. To ensure a smoother installation and avoid conflicting dependencies, it's highly recommended to create a dedicated environment for PPanGGOLiN:

```bash
# Install into a new conda environment
conda create -n ppanggo -c conda-forge -c bioconda ppanggolin

# Check PPanGGOLiN install
conda activate ppanggo
ppanggolin --version
```

```{tip}
  For quicker installation and potential resolution of conflicting dependencies, consider using Mamba, an efficient alternative to conda.
```

```{note}
Supported python version are 3.8, 3.9 and 3.10
```


## Installing from Source Code (GitHub)

### Within a Conda Environment

A straightforward method to install PPanGGOLiN from the source code is by utilizing a conda environment that includes all the necessary dependencies.

**1. Clone the PPanGGOLiN Repository**

```bash
git clone https://github.com/labgem/PPanGGOLiN.git
cd PPanGGOLiN
```

**2. Installing Dependencies with a Conda Environment File**

Install PPanGGOLiN dependencies listed in the [ppanggolin_env.yaml](../../ppanggolin_env.yaml) file, located at the root of the repository, using conda:

```bash
conda create -n ppanggolin_source -f ppanggolin_env.yaml
```

**3. Installing PPanGGOLiN**

Finally, install PPanGGOLiN using **pip**:

```bash
pip install .
```

### Manual Installation

If you prefer to install PPanGGOLiN without using conda, follow these steps:

```{warning}
Please note that this method involves manual installation. While assistance is available, complete troubleshooting support may not be provided. We strongly recommend using the conda installation method whenever possible.
```

**1. Clone the PPanGGOLiN Repository**

```bash
git clone https://github.com/labgem/PPanGGOLiN.git
cd PPanGGOLiN
```

**2. Installing PPanGGOLiN Dependencies**

To ensure the tool functions correctly, you need to install all dependencies listed in the [ppanggolin_env.yaml](../../ppanggolin_env.yaml) file.

**2.1. Installing Python Dependencies**

Create a file named `requirements.txt` and add the following contents:

```text
tqdm>=4.64
pytables>=3.7
pyrodigal>=3.0.1
networkx>=3.0
dataclasses>=0.8
scipy>=1.10.0
plotly>=4.14.3
gmpy2>=2.1.2
pandas>=2.0
colorlover>=0.3
numpy>=1.24
bokeh>=2.4.2,<3
```

Then, use **pip** to install these dependencies:

```bash
python3 -m pip install -r requirements.txt
```

```{warning}
Ensure you are using Python version 3.8 or higher.
```

**2.2. Installing Required Software**

Install the following non-Python software:

- [MMSeqs2>=13.45111](https://github.com/soedinglab/MMseqs2/wiki#installation)
- [Aragorn>=1.2.41](http://www.ansikte.se/ARAGORN/Downloads/)
- [Infernal>=1.1.4](http://eddylab.org/infernal/)
- [MAFFT>=7.505](https://mafft.cbrc.jp/alignment/software/)

```{note}
- MMSeqs2 is crucial for gene clustering, while Aragorn and Infernal are used for genome annotation.
- MAFFT is utilized in the `ppanggolin msa` command for multiple sequence alignment.
- Skip installing Aragorn, Infernal, or MAFFT if you do not require their specific features.
```

**3. Installing PPanGGOLiN**

Finally, install PPanGGOLiN using **pip**:

```bash
pip install .
```

## Development Version

If you wish to utilize the development version of PPanGGOLiN, you can access the 'dev' branch on GitHub. Please note that while its functionality is not guaranteed, it typically works most of the time.

Follow these steps to obtain and install the development version:

**1. Clone the Repository**

Clone the 'dev' branch of the repository onto your local machine:

```bash
git clone --branch dev https://github.com/labgem/PPanGGOLiN.git
cd PPanGGOLiN
```

**2. Install Dependencies**

Ensure you have all the necessary dependencies installed. Refer to the [installation instructions above](#installing-from-source-code-github) for guidance on installing dependencies.

**3. Install PPanGGOLiN**

Once dependencies are installed, proceed to install PPanGGOLiN using **pip**:

```bash
pip install .
```

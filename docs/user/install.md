# Installation

```{note}
Supported Python versions are 3.9, 3.10, 3.11 and 3.12
```

## Installing PPanGGOLiN with Conda (Recommended)

The recommended way to install PPanGGOLiN is via Conda from the Bioconda channel. To ensure smooth installation and avoid dependency conflicts, it's best to create a dedicated Conda environment for PPanGGOLiN:

```bash
# Create a new Conda environment and install PPanGGOLiN
conda create -n ppanggo -c conda-forge -c bioconda ppanggolin

# Activate the environment
conda activate ppanggo

# Verify the installation
ppanggolin --version
```

## Installing PPanGGOLiN from PyPI

Starting from version 2.2.2, PPanGGOLiN is available on PyPI and can be installed using pip:

```bash
pip install ppanggolin

# Verify the installation
ppanggolin --version
```

```{warning}
While this method installs PPanGGOLiN and its Python dependencies, it does not install the non-Python dependencies required for many commands within PPanGGOLiN. Be sure to install them separately to ensure full functionality. See the [Installing PPanGGOLiN Dependencies](#install_dependencies) section for more details.
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
conda env create -n ppanggolin_source -f ppanggolin_env.yaml
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

**2. Installing PPanGGOLiN Dependencies** <a name="install_dependencies"></a>

To ensure the tool functions correctly, you need to install all dependencies listed in the [ppanggolin_env.yaml](../../ppanggolin_env.yaml) file.

Next, install the following non-Python software:

- [MMSeqs2>=13.45111](https://github.com/soedinglab/MMseqs2/wiki#installation)
- [Aragorn>=1.2.41](http://www.ansikte.se/ARAGORN/Downloads/)
- [Infernal>=1.1.4](http://eddylab.org/infernal/)
- [MAFFT>=7.505](https://mafft.cbrc.jp/alignment/software/)

```{note}
- MMSeqs2 is crucial for gene clustering, while Aragorn and Infernal are used for genome annotation.
- MAFFT is utilized in the `ppanggolin msa` command for multiple sequence alignment.
- Skip installing Aragorn, Infernal, or MAFFT if you do not require their specific features.
```

For Debian, Ubuntu, Linux Mint, Kali Linux, and all other distributions using apt you can run the following command:

```bash
sudo apt install mmseqs2 infernal aragorn mafft
```


**3. Installing PPanGGOLiN and its Python Dependencies**

PPanGGOLiN's Python dependencies are specified in the `pyproject.toml` file under the optional dependencies category named `python_deps`. This configuration file is situated at the root of the PPanGGOLiN repository.

To install PPanGGOLiN along with its Python dependencies, you can use the following pip command:

```bash
pip install .
```


## Development Version

If you wish to utilize the development version of PPanGGOLiN, you can access the 'dev' branch on GitHub. Please note that while its functionality is not guaranteed, it typically works most of the time.

Follow these steps to obtain and install the development version:

**1. Clone the Repository**

Clone the `dev` branch of the repository onto your local machine:

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

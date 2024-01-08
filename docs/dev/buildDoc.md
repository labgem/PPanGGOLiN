# Building the Documentation

This section provides guidelines for building the PPanGGOLiN documentation locally.

## Setting Up the Environment

Before proceeding, ensure that you have installed PPanGGOLiN from the source code. For detailed instructions, refer to [this section](../user/install.md#installing-from-source-code-github).

The necessary packages to build the documentation are listed in the 'requirements.txt' file located in the `doc/` folder.

```bash
pip install -r docs/requirements.txt
```

Alternatively, the same list of requirements is available in the [pyproject.toml file](../../pyproject.toml), which allows for automatic installation using `pip`.

```shell
# Replace '/path/to/ppanggolin/' with your actual path
pip install /path/to/ppanggolin/[doc]
```

## Building Documentation with Sphinx
### Build and produce an html

Building the documentation is as simple as :


```bash
# Replace '/path/to/ppanggolin/' with your actual path
cd /path/to/ppanggolin/docs/
sphinx-build -b html . build/
```
You can also use the makefile as follow


```bash
# Replace '/path/to/ppanggolin/' with your actual path
cd /path/to/ppanggolin/docs/
make html
```


### Build with autobuild 

You can visualize your modifications in real-time using **sphinx-autobuild**, a tool previously installed.

```shell
cd $PPanGGOLiN/docs
sphinx-autobuild . build/
# Copy the server address, for example: http://127.0.0.1:8000
# Paste the address in your browser
```

```{note}
The package [readthedocs-sphinx-search](https://readthedocs-sphinx-search.readthedocs.io/en/latest/) enables "search as you type" functionality for docs hosted on Read the Docs. Please note that it only functions on the ReadTheDocs website. `[INFO] Docs are not being served on Read the Docs, readthedocs-sphinx-search will not work.`
```

### Editing or Adding Documentation


To modify existing documentation:

1. **Navigate to the Document**: Go to the file you wish to edit and make necessary changes.

To add a new page:

1. **Create Markdown File**: Place the new markdown file in the relevant folder within the 'docs' directory—'user' for user documentation or 'dev' for developer documentation.
2. **Update Table of Contents (TOC)**: Add a reference to the newly added file in the 'index.md' file at the root of the docs folder under the 'user' or 'dev' TOC tree.



(add-api-doc)=
#### Update API documentation
The API reference documentation is automatically updated as it is build each time the doc is build by sphinx.
In case you added a new file in the ppanggolin code base 

To update the API documentation and keep the automatic update when a new package, module, submodules is added follow the
next lines:
```shell
sphinx-apidoc -o api $PPanGGOLiN/ppanggolin -f
```
```{note}
*sphinx-apidoc* will generate ReStructeredText files. You need to convert them in markdown. For this follow the guides 
[here](#rst2md)
```


### Creating a New Documentation from Scratch

This section documents how the current documentation has been created. 

#### Quickstart with Sphinx

To start the documentation process from scratch, follow these steps to either rename the existing documentation or provide a new name for the upcoming documentation.

```shell
DOCS=path/to/PPanGGOLiN/docs
sphinx-quickstart $DOCS
```

Upon executing the command, you will be prompted with a series of settings in order to setup the new documentation folder.

We used so far the default settings as follow:

- Separate source and build directories (y/n) [n]: **n**

-  Project name: **PPanGGOLiN**
-  Author name(s): **Your name**
-  Project release []: **The current version of PPanGGOLiN**



#### Configuration file

Locate the conf.py file within the docs directory. You can modify this file similarly to the adjustments made in the current `conf.py` file.

(rst2md)=
#### ReStructeredText to markdown

reStructuredText (rst) is the default plaintext markup language used by both Docutils and Sphinx. 
Despite being more comprehensive, it's considered slightly older and less user-friendly compared to Markdown.

We have decided to use Markdown (md) instead of reStructuredText for our documentation
We will use [MyST](https://mystmd.org/guide) to translate RST files to Markdown while preserving all features provided by reStructuredText.

For this we will need to install the package `rst-to-myst`.


```shell
pip install rst-to-myst

rst2myst convert index.rst

# remove rst file(s)
rm index.rst 
```

#### User documentation
Here are some general guidelines to write user documentation:

1. **Topic/Command Separation**: Create individual files for each topic or command, offering explicit explanations of the feature's functionality.
    - **Enhance with Examples**: Include example code snippets and output figures or initial lines of output files wherever applicable.

2. **Clarity and Precision**: Strive for utmost clarity by defining acronyms and jargon used within the documentation.


#### API documentation

To generate the API documentation using the docstrings in your code, follow these steps:

1. Using `sphinx-apidoc`:

Generate the API documentation files with the `sphinx-apidoc` command:

```bash
# Generate API doc files
sphinx-apidoc -o api $PPanGGOLiN/ppanggolin
```

This command creates an 'api' folder containing the skeleton of the API reference pages.

2. Translating to MyST with `rst-to-myst`:

Translate these generated files into MyST markdown using `rst-to-myst`:

```bash
# Translate them into MyST
rst2myst convert api/*.rst

# Remove remaining RST files
rm api/*.rst
```

```{tip}
With the "sphinx.ext.autosectionlabel", you will certainly get multiple warning for duplicate label. 
To remove them you have to remove or modify the label in one of the cited file.
```
```{tip}
When you use "sphinx-apidoc" a modules.md file is created but he is not used. we advice to removed it to prevent warning.
```
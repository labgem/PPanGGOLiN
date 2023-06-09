# Using conda

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

# Using github

If you want to use the development version you can use the 'master' branch on Github. While it is not guaranteed to work, it should most of the time.

You need to clone the repository on your computer, then you will need to have all of the dependencies (they are listed in the [requirements](https://github.com/labgem/PPanGGOLiN/blob/master/requirements.txt) file at the root of the repositories). You can install them using conda as such : 
` conda install --file requirements.txt `

Then you will need to install ppanggolin. You can use pip while in your conda environment which is the preferred way of doing. For this, place yourself at the root of the repository and run : 

`pip install .`

It should install PPanGGOLiN in your conda environment.
Alternatively, you can use the setup.py file directly, as such  : 

`python setup.py install .`
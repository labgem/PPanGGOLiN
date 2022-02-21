#!/usr/bin/env python3

import setuptools
import os

from distutils.extension import Extension

NEM_DIR_PATH = "ppanggolin/nem/NEM/"

if __name__ == "__main__":
    setuptools.setup(
        name="ppanggolin",
        version=open(os.path.join(os.path.dirname(__file__), "VERSION")).read().rstrip(),
        url="https://github.com/labgem/PPanGGOLiN",
        description="Pangenome analysis suite",
        packages=setuptools.find_packages(),
        setup_requires=["cython"],
        install_requires=[],
        package_data={'': ['rRNA_DB/*cm*']},
        classifiers=["Environment :: Console",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
                     "Natural Language :: English",
                     "Operating System :: POSIX :: Linux",
                     "Programming Language :: Python :: 3",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"],
        entry_points={"console_scripts": ["ppanggolin = ppanggolin.main:main"]},
        ext_modules=[Extension(
            extra_compile_args=['-fcommon'],
            name="nem_stats",
            sources=[NEM_DIR_PATH + 'nem_stats.pyx',
                     NEM_DIR_PATH + 'nem_exe.c',
                     NEM_DIR_PATH + 'nem_alg.c',
                     NEM_DIR_PATH + 'nem_nei.c',
                     NEM_DIR_PATH + 'nem_mod.c',
                     NEM_DIR_PATH + 'nem_rnd.c',
                     NEM_DIR_PATH + 'lib_io.c',
                     NEM_DIR_PATH + 'nem_hlp.c',
                     NEM_DIR_PATH + 'genmemo.c'],
            include_dirs=[NEM_DIR_PATH])]
    )

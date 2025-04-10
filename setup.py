#!/usr/bin/env python3

from setuptools import Extension, setup

NEM_DIR_PATH = "ppanggolin/nem/NEM/"

setup(
    ext_modules=[
        Extension(
            # Force c99, NEM uses features deprecated in recent C standards
            extra_compile_args=["-fcommon", "-Wno-int-conversion", "-std=c99"],
            name="nem_stats",
            sources=[
                NEM_DIR_PATH + "nem_stats.pyx",
                NEM_DIR_PATH + "nem_exe.c",
                NEM_DIR_PATH + "nem_alg.c",
                NEM_DIR_PATH + "nem_nei.c",
                NEM_DIR_PATH + "nem_mod.c",
                NEM_DIR_PATH + "nem_rnd.c",
                NEM_DIR_PATH + "lib_io.c",
                NEM_DIR_PATH + "nem_hlp.c",
                NEM_DIR_PATH + "genmemo.c",
            ],
            include_dirs=[NEM_DIR_PATH],
        )
    ]
)

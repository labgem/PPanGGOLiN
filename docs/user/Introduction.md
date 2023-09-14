# Introduction

PPanGGOLiN is built as a software suite to analyse groups of genomes potentially belonging to the same taxonomic group (subspecies, species, genus...). 
It is an open source CeCiLL-licensed software, implemented in C and python3 and works under Linux and MacOS.
It is built to use multiple cores and is able to handle tens of thousands of genomes in a single analysis with (relatively) reasonable CPU and RAM usage.


The software has a minimal one-liner that can be used to generate most of the outputs of PPanGGOLiN without parameter tuning.
While it should work fine at the species or subspecies level, it will yield indigent results with more divergent genomes, and working with PPanGGOLiN on such datasets will require parameter tuning.


If you want to tune your analysis, the PPanGGOLiN suite can be used through multiple subcommands that allow tuning parameters for each step of the analysis.
The results of each analysis are saved in the same file, which is in the HDF-5 file format (with the .h5 extension). 
When you run an analysis using this file as input, the results of that analysis will be added to the file to supplement the data that are already stored in it. 
The idea behind this is that you can store and manipulate your pangenomes with PPanGGOLiN by using this file only. It will keep all the information about what was done, all the parameters used, and all the pangenome's content.


Keeping it this way makes it possible to use the file as a pangenome reference database for a taxonomic group for some applications such as evaluating MAG completeness or contamination, or comparing a genome with a pangenome to assign partitions to the genome's genes. Other applications are currently in development.
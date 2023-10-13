
Once partitions have been computed, you can predict the regions of genome plasticity. 
This subcommand's options are about tuning parameters for the analysis. Details about each parameter can be found in the related [article](https://doi.org/10.1093/bioinformatics/btaa792).

You can do it as such:

`ppanggolin rgp -p pangenome.h5`

This will predict RGPs and store results in the HDF5 file. If you want a list of RGPs for each genome, you can use `ppanggolin write -p pangenome.h5 --regions --output MYOUTPUTDIR`. It will provide the file 'plastic regions' whose format is described [here](https://github.com/labgem/PPanGGOLiN/wiki/Outputs#plastic-regions)
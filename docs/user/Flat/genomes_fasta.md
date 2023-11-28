<!-- ### Adding Fasta Sequences into GFF and proksee JSON map Files -->

PPanGGOLiN allows the incorporation of fasta sequences into GFF files and proksee JSON map files. This integration with Proksee provides access to various tools that rely on DNA sequences, including the construction of GC% and GC skew profiles, and conducting blast searches for example.


Since PPanGGOLiN does not retain genomic sequences, it is necessary to provide the original genomic files used to construct the pangenome through either the `--anno` or `--fasta` argument. These arguments mirror those used in workflow commands (`workflow`, `all`, `panrgp`, `panmodule`) and the `annotate` command.

- `--anno`: This option requires a tab-separated file containing organism names and the corresponding GFF/GBFF filepaths of their annotations. If `--anno` is utilized, GFF files should include fasta sequences.

- `--fasta`: Use this option with a tab-separated file that lists organism names alongside the filepaths of their genomic sequences in fasta format.


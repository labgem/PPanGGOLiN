# Multiple Sequence Alignment

This command is available from 1.1.103 and on.
It is used to call [mafft](https://mafft.cbrc.jp/alignment/software/) with default options to compute MSA of any partition of the pangenome. Using multiple cpus is recommended as it is quite demanding in computational resources.

By default it will write the strict 'core' (genes that are present in absolutely all genomes) and remove any duplicated genes. Beware however that, if you have many genomes (over 1000), the core will likely be either very small or even empty if you have fragmented genomes.

It will write one MSA for each family. You can then provide the directory where the MSA are written to [IQ-TREE](https://github.com/Cibiv/IQ-TREE) for example, to do phylogenetic analysis.

### partitions

You can change the partition which is written, by using the --partition option.
`ppanggolin msa -p pangenome.h5 --partition persistent` for example will compute MSA for all the persistent gene families.

Supported partitions are core, persistent, shell, cloud, softcore, accessory. If you wish to have additional filters, you can raise an issue with your demand, or write a PR directly, most possibilites should be quite straightforward to add.

### source

You can specify whether to use dna or protein sequences for the MSA by using --source. It uses protein sequences by default.

`ppanggolin msa -p pangenome.h5 --source dna`

### phylo

It is also possible to write a single whole genome MSA file, which many phylogenetic softwares accept as input, by using the --phylo option as such:

`ppanggolin msa -p pangenome.h5 --phylo`

This will contatenate all of the family MSA into a single MSA, with one sequence for each genome.
# Align external genes to a pangenome


`ppanggolin align` is a command made to use a pangenome as a reference to get information about a set of sequences of interest. It requires a HDF-5 file of a previously computed pangenome as input, as well as either a set of sequences, nucleotides or proteins, as a .fasta file.
The command will use MMseqs to compare the given input sequences to the pangenome gene family representatives, and will assign a gene family to each input sequence, if there is one which is close enough (the 'closeness' can be defined by some of the command parameters). If there are multiple assignable families, the closest one in bitscore is chosen.

The default way of using this command is the following:

`ppanggolin align -p pangenome.h5 -o MYOUTPUTDIR --sequences MY_SEQUENCSE_OF_INTEREST.fasta`


By default the command if successful will always write two files:

- The file 'sequences_partition_projection.tsv' which is a .tsv file with two columns which indicates the partition of the most similar gene family of the pangenome to which the given input sequence is closest to. It follows the following format:

| column | description |
|--------|-------------|
| input | the header of the sequence in the given .fasta file|
|partition| predicted partition based on the most similar gene family, or 'cloud' if there are no similar enough gene family|

- The file 'input_to_pangenome_associations.blast-tab' is a .tsv file which follows the tabular blast format which many alignment softwares (such as blast, diamond, mmseqs etc.) use, with two additional columns: the length of query sequence which was aligned, and the length of the subject sequence which was aligned (provided with qlen and slen with the softwares I previously named). You can find a detailed description of the format in [this blog post](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) for example (and there are many other descriptions of this format on internet, if you search for 'tabular blast format'). The query are the provided sequences, and the subjet are the pangenome gene families.

As options, you can also write additional files which give you alternative informations. If your pangenome has predicted RGP and spots (see [Regions of Genome Plasticity](RGP/rgpAnalyses.md) if you do not know what those are) 
you can use `--getinfo` as such:

`ppanggolin align -p pangenome.h5 -o MYOUTPUTDIR --sequences MY_SEQUENCSE_OF_INTEREST.fasta --getinfo`

`--getinfo` will list known spots and RGPs where the gene families similar to your proteins of interest are found. They will be listed if they are in the RGPs themselves OR if they are bordering it (that is, if they are within 3 persistent genes of the RGP).
The written file will be called 'info_input_seq.tsv', and follows the following format:

| column | description |
|--------|-------------|
| input | the header of the sequence in the given .fasta file|
|family| the id of the family the input sequence was assigned to|
|partition| predicted partition based on the most similar gene family, or 'cloud' if there are no similar enough gene family|
|spot_list_as_member| the list of spots in which the sequence is found, as a member of the spot (it is included in it)|
|spot_list_as_border| the list of spots in which the sequence is found as a bordering gene|
|rgp_list| the list of RGP in which the sequence is found|

You can use `--draw_related` as such:
`ppanggolin align -p pangenome.h5 -o MYOUTPUTDIR --sequences MY_SEQUENCSE_OF_INTEREST.fasta --draw_related`

Which will draw all of the spots where the gene families similar to your proteins of interest are found, writing 3 files, one figure, one .gexf file and one .tsv file. This option is basically using what is described in the [`draw --spots`](RGP/rgpOutputs.md#draw-spots) part of the documentation.

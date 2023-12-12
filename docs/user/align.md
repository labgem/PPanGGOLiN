# Align external genes to a pangenome


The PPanGGOLiN `align` command allows to use a pangenome as a reference to get information about a set of sequences of interest. It requires a previously computed pangenome in HDF-5 format as input, along with a `.fasta` file containing either nucleotide or protein sequences.
The command utilizes MMseqs to compare input sequences to representatives of the pangenome gene family. It assigns a gene family to each input sequence if there is one that is sufficiently similar (as defined by command parameters). If multiple families are assignable, the one with the highest bitscore is selected.

This command is used as follows:
```bash
ppanggolin align -p pangenome.h5 -o MYOUTPUTDIR --sequences MY_SEQUENCSE_OF_INTEREST.fasta
```

## Output files
By default the command creates two output files:

### 1. 'sequences_partition_projection.tsv'

'sequences_partition_projection.tsv' is a .tsv file with two columns that indicates the partition of the most similar gene family in the pangenome to which the given input sequence is closest. It follows the following format:

| column | description |    
|--------|-------------|
| input | the header of the sequence in the given .fasta file|
|partition| predicted partition based on the most similar gene family, or 'cloud' if there are <br> no similar enough gene family|


### 2. 'input_to_pangenome_associations.blast-tab'

'input_to_pangenome_associations.blast-tab' is a .tsv file that follows the tabular blast format which many alignment softwares (such as blast, diamond, mmseqs etc.) use, with two additional columns: the length of query sequence which was aligned, and the length of the subject sequence which was aligned (provided with qlen and slen with the softwares I previously named). You can find a detailed description of the format in [this blog post](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) for example (and there are many other descriptions of this format on internet, if you search for 'tabular blast format'). The query are the provided sequences, and the subjet are the pangenome gene families.


### 3. Optional outputs 

Optionally, you can also write additional files that provide alternative information. If RGP and spots have been predicted in your pangenome (see [Regions of Genome Plasticity](RGP/rgpAnalyses.md) if you do not know what those are) 
you can use `--getinfo` as such:

```bash
ppanggolin align -p pangenome.h5 -o MYOUTPUTDIR --sequences MY_SEQUENCSE_OF_INTEREST.fasta --getinfo
```


`--getinfo` will list known spots and RGPs where the gene families similar to your proteins of interest are found. They will be listed if they are in the RGPs themselves OR if they are bordering it (that is, if they are within 3 persistent genes of the RGP).
The written file will be called 'info_input_seq.tsv', and follows the following format:

| column | description |
|--------|-------------|
| input | the header of the sequence in the given .fasta file|
|family| the id of the family the input sequence was assigned to|
|partition| predicted partition based on the most similar gene family, or 'cloud' if <br> there are no similar enough gene family|
|spot_list_as_member| the list of spots in which the sequence is found, as a member of the spot <br> (it is included in it)|
|spot_list_as_border| the list of spots in which the sequence is found as a bordering gene|
|rgp_list| the list of RGP in which the sequence is found|

You can use `--draw_related` as such:
```bash
ppanggolin align -p pangenome.h5 -o MYOUTPUTDIR --sequences MY_SEQUENCSE_OF_INTEREST.fasta --draw_related
```


It will draw all of the spots where the gene families similar to your proteins of interest are found, writing 3 files, one figure, one .gexf file and one .tsv file. This option is basically using what is described in the [`draw --spots`](RGP/rgpOutputs.md#draw-spots) part of the documentation.

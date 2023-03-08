#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import tempfile
import subprocess
import time
from multiprocessing import get_context

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir, restricted_float
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.genetic_codes import genetic_codes


def is_single_copy(fam, dup_margin):
    """
    Check if a gene family can be considered 'single copy' or not
    
    :param fam: GeneFamily object
    :param dup_margin: maximal number of genomes in which the gene family can have multiple members and still be considered a 'single copy' gene family
    """
    nb_multi = 0
    for gene_list in fam.get_org_dict().values():
        if len(gene_list) > 1:
            nb_multi += 1
    dup_ratio = nb_multi / len(fam.organisms)
    if dup_ratio < dup_margin:
        return True
    return False


def getFamiliesToWrite(pangenome, partition_filter, soft_core=0.95, dup_margin=0.95, single_copy=True):

    """
    Get families corresponding to the given partition

    :param pangenome: Partitioned pangenome
    :param partition_filter: choice of partition to compute Multiple Sequence Alignement of the gene families
    :param soft_core: Soft core threshold to use
    :param dup_margin: maximal number of genomes in which the gene family can have multiple members and still be considered a 'single copy' gene family
    :param single_copy: Use "single copy" (defined by dup_margin) gene families only

    :return: set of families unique to one partition
    """
    fams = set()
    nb_org = pangenome.number_of_organisms()

    if partition_filter == "all":
        return set(pangenome.gene_families)
    if partition_filter in ["persistent", "shell", "cloud"]:
        for fam in pangenome.gene_families:
            if fam.named_partition == partition_filter:
                if single_copy and is_single_copy(fam, dup_margin):
                    fams.add(fam)
                elif not single_copy:
                    fams.add(fam)
    elif partition_filter in ["core", "accessory", "softcore"]:
        if partition_filter == "core":
            for fam in pangenome.gene_families:
                if len(fam.organisms) == nb_org:
                    if single_copy and is_single_copy(fam, dup_margin):
                        fams.add(fam)
                    elif not single_copy:
                        fams.add(fam)
        elif partition_filter == "accessory":
            for fam in pangenome.gene_families:
                if len(fam.organisms) < nb_org:
                    if single_copy and is_single_copy(fam, dup_margin):
                        fams.add(fam)
                    elif not single_copy:
                        fams.add(fam)
        elif partition_filter == "softcore":
            for fam in pangenome.gene_families:
                if len(fam.organisms) >= nb_org * soft_core:
                    if single_copy and is_single_copy(fam, dup_margin):
                        fams.add(fam)
                    elif not single_copy:
                        fams.add(fam)
    return fams


def translate(seq: str, code: dict):
    """translates the given dna sequence with the given translation table

    :param seq: given dna sequence
    :param code: translation table corresponding to genetic code to use

    :return: protein sequence
    """
    # code:  https://www.bioinformatics.org/sms/iupac.html
    start_table = code["start_table"]
    table = code["trans_table"]

    if len(seq) % 3 == 0:
        protein = start_table[seq[0: 3]]
        for i in range(3, len(seq), 3):
            codon = seq[i: i + 3]
            try:
                protein += table[codon]
            except KeyError:  # codon was not planned for. Probably can't determine it.
                protein += 'X'  # X is for unknown
    else:
        raise IndexError("Given sequence length modulo 3 was different than 0, which is unexpected.")
    return protein


def write_fasta_families(family: GeneFamily, tmpdir: tempfile.TemporaryDirectory, code_table: dict,
                         source: str = 'protein', use_gene_id: bool = False):
    """Write fasta files for each gene family

    :param family: gene family to write
    :param tmpdir: path to temporary directory
    :param source: indicates whether to use protein or dna sequences to compute the msa
    :param use_gene_id: Use gene identifiers rather than organism names for sequences in the family MSA
    :param code_table: Genetic code to use

    :return: path to fasta file
    """
    # have a directory for each gene family, to make deletion of tmp files simpler

    f_name = tmpdir.name + "/" + family.name + ".fasta"
    f_obj = open(f_name, "w")
    # get genes that are present in only one copy for our family in each organism.
    single_copy_genes = []
    for _, genes in family.get_org_dict().items():
        if len(genes) == 1:
            single_copy_genes.extend(genes)

    for gene in single_copy_genes:
        if use_gene_id:
            f_obj.write('>' + gene.ID + "\n")
        else:
            f_obj.write('>' + gene.organism.name + "\n")
        if source == "dna":
            f_obj.write(gene.dna + '\n')
        elif source == "protein":
            f_obj.write(translate(gene.dna, code_table) + "\n")
        else:
            raise Exception("Unknown sequence source given (expected 'dna' or 'protein')")
    f_obj.flush()

    return f_name


def launch_mafft(fname, output, fam_name):
    """
    Compute the MSA with mafft

    :param fname: family gene sequence in fasta
    :param output: directory to save alignment
    :param fam_name: Name of the gene family
    """
    outname = output + "/" + fam_name + ".aln"
    cmd = ["mafft", "--thread", "1", fname]
    logging.getLogger().debug("command: " + " ".join(cmd))
    subprocess.run(cmd, stdout=open(outname, "w"), stderr=subprocess.DEVNULL, check=True)


def launch_multi_mafft(args):
    """ Allow to launch mafft in multiprocessing

    :param args: Pack of argument for launch_mafft

    :return: Organism object for pangenome
    """
    launch_mafft(*args)


def compute_msa(families: set, output: str, tmpdir: str, cpu: int = 1, source: str = "protein",
                use_gene_id: bool = False, code: int = 11, disable_bar: bool = False):
    """
    Compute MSA between pangenome gene families

    :param families: Set of families specific to given partition
    :param output: output directory name for families alignment
    :param cpu: number of available core
    :param tmpdir: path to temporary directory
    :param source: indicates whether to use protein or dna sequences to compute the msa
    :param use_gene_id: Use gene identifiers rather than organism names for sequences in the family MSA
    :param code: Genetic code to use
    :param disable_bar: Disable progress bar
    """
    newtmpdir = tempfile.TemporaryDirectory(dir=tmpdir)

    write_total = 0
    args = []
    logging.getLogger().info("Preparing input files for MSA...")
    code_table = genetic_codes(str(code))

    for family in tqdm(families, unit="family", disable=disable_bar):
        start_write = time.time()
        fname = write_fasta_families(family, newtmpdir, code_table, source, use_gene_id)
        write_total = write_total + (time.time() - start_write)
        args.append((fname, output, family.name))

    logging.getLogger().info("Computing the MSA ...")
    bar = tqdm(range(len(families)), unit="family", disable=disable_bar)
    with get_context('fork').Pool(cpu) as p:
        for _ in p.imap_unordered(launch_multi_mafft, args):
            bar.update()
    bar.close()


def write_whole_genome_msa(pangenome: Pangenome, families: set, phylo_name: str, outname: str,
                           use_gene_id: bool = False):
    """
    Writes a whole genome msa file for additional phylogenetic analysis

    :param pangenome: Pangenome object
    :param families: Set of families specific to given partition
    :param phylo_name: output file name for phylo alignment
    :param outname: output directory name for families alignment
    :param use_gene_id: Use gene identifiers rather than organism names for sequences in the family MSA
    """
    phylo_dict = {}
    for org in pangenome.organisms:
        phylo_dict[org.name] = ""
    for fam in families:
        missing_genomes = set(phylo_dict.keys())
        fin = open(outname + "/" + fam.name + ".aln", "r")
        genome_id = ""
        seq = ""
        curr_len = 0
        dup_gene = 0
        curr_phylo_dict = {}

        for line in fin:
            if line.startswith('>'):
                if genome_id != "":
                    if genome_id not in missing_genomes:
                        dup_gene += 1
                        # duplicated genes. Replacing them with gaps.
                        curr_phylo_dict[genome_id] = "-" * curr_len
                    else:
                        curr_phylo_dict[genome_id] = seq
                        missing_genomes -= {genome_id}
                        curr_len = len(seq)
                if use_gene_id:
                    genome_id = pangenome.get_gene(line[1:].strip()).organism.name
                else:
                    genome_id = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
        if genome_id != "":
            if genome_id not in missing_genomes:
                # duplicated genes. Replacing them with gaps.
                curr_phylo_dict[genome_id] = "-" * curr_len
            else:
                curr_phylo_dict[genome_id] = seq
                curr_len = len(seq)
        fin.close()
        for genome in missing_genomes:
            curr_phylo_dict[genome] = "-" * curr_len

        for key, val in curr_phylo_dict.items():
            phylo_dict[key] += val

    fout = open(phylo_name, "w")
    for key, val in phylo_dict.items():
        fout.write(">" + key + "\n")
        fout.write(val + "\n")
    fout.close()


def writeMSAFiles(pangenome, output, cpu=1, partition="core", tmpdir="/tmp", source="protein", soft_core=0.95,
                  phylo=False, use_gene_id=False, translation_table="11", dup_margin = 0.95, single_copy=True, force=False, disable_bar=False):

    """
    Main function to write MSA files

    :param pangenome: Pangenome object with partition
    :param output: Path to output directory
    :param cpu: number of available core
    :param partition: choice of partition to compute Multiple Sequence Alignement of the gene families
    :param tmpdir: path to temporary directory
    :param source: indicates whether to use protein or dna sequences to compute the msa
    :param soft_core: Soft core threshold to use
    :param phylo: Writes a whole genome msa file for additional phylogenetic analysis
    :param use_gene_id: Use gene identifiers rather than organism names for sequences in the family MSA
    :param translation_table: Translation table (genetic code) to use.
    :param dup_margin: maximal number of genomes in which the gene family can have multiple members and still be considered a 'single copy' gene family
    :param single_copy: Use "single copy" (defined by dup_margin) gene families only
    :param force: force to write in the directory
    :param disable_bar: Disable progress bar
    """
    need_partitions = False
    if partition in ["persistent", "shell", "cloud"]:
        need_partitions = True

    outname = output + f"/msa_{partition}_{source}/"
    mk_outdir(outname, force=force)

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_partitions=need_partitions,
                         need_gene_sequences=True, disable_bar=disable_bar)
    logging.getLogger().info(f"Doing MSA for {partition} families...")
    families = getFamiliesToWrite(pangenome, partition_filter=partition, soft_core=soft_core, dup_margin=dup_margin, single_copy=single_copy)

    # check that the code is similar than the one used previously, if there is one
    if 'translation_table' in pangenome.parameters["cluster"]:
        if pangenome.parameters["cluster"]["translation_table"] != translation_table:
            logging.getLogger().warning("The translation table used during clustering "
                                        f"('{pangenome.parameters['cluster']['translation_table']}') "
                                        f"is different than the one provided now ('{translation_table}')")
    code = translation_table

    compute_msa(families, outname, cpu=cpu, tmpdir=tmpdir, source=source, use_gene_id=use_gene_id, code=code,
                disable_bar=disable_bar)
    logging.getLogger().info(f"Done writing all {partition} MSA in: {outname}")

    if phylo:
        logging.getLogger().info("Writing the whole genome msa file")
        if partition == "softcore":
            phylo_name = output + f"/{partition}_{soft_core}_genome_alignment.aln"
        else:
            phylo_name = output + f"/{partition}_genome_alignment.aln"
        write_whole_genome_msa(pangenome, families, phylo_name, outname, use_gene_id=use_gene_id)
        logging.getLogger().info(f"Done writing the {partition} genome alignment in: '{phylo_name}'")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    writeMSAFiles(pangenome, args.output, cpu=args.cpu, partition=args.partition, tmpdir=args.tmpdir,
                  source=args.source, soft_core=args.soft_core, phylo=args.phylo, use_gene_id=args.use_gene_id,
                  translation_table=args.translation_table, dup_margin=args.dup_margin,
                  single_copy=args.single_copy, force=args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("msa", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_msa(parser)
    return parser


def parser_msa(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of msa command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="The following arguments are required :")
    required.add_argument('-p', '--pangenome', required=False, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title="Optional arguments. Indicating 'all' writes all elements. "
                                               "Writing a partition ('persistent', 'shell', 'cloud', 'core' or "
                                               "'accessory') write the elements associated to said partition.")
    # could make choice to allow customization
    optional.add_argument("--soft_core", required=False, type=restricted_float, default=0.95,
                          help="Soft core threshold to use if 'softcore' partition is chosen")
    optional.add_argument("--dup_margin", required=False, type=restricted_float, default=0.05,
                          help="minimum ratio of organisms in which the family must have multiple genes "
                               "for it to be considered 'duplicated'")
    optional.add_argument("--single_copy", required=False, action="store_true", default=False,
                          help="Use report gene families that are considered 'single copy', for details see "
                               "option --dup_margin")
    optional.add_argument("--partition", required=False, default="core",
                          choices=["all", "persistent", "shell", "cloud", "core", "accessory", 'softcore'],
                          help="compute Multiple Sequence Alignement of the gene families in the given partition")
    optional.add_argument("--source", required=False, default="protein", choices=["dna", "protein"],
                          help="indicates whether to use protein or dna sequences to compute the msa")
    optional.add_argument("--phylo", required=False, action='store_true',
                          help="Writes a whole genome msa file for additional phylogenetic analysis")
    optional.add_argument("--use_gene_id", required=False, action='store_true',
                          help="Use gene identifiers rather than organism names for sequences in the family MSA"
                               " (organism names are used by default)")
    optional.add_argument("--translation_table", required=False, default=11, type=int,
                          help="Translation table (genetic code) to use.")

    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_msa(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(),
                        help="directory for storing temporary files")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())

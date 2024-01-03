#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import tempfile
import subprocess
import time
from multiprocessing import get_context
from pathlib import Path
from typing import Dict, Set, List, Tuple

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import mk_outdir, restricted_float
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.genetic_codes import genetic_codes


def is_single_copy(family: GeneFamily, dup_margin: float = 0.95) -> bool:
    """
    Check if a gene family can be considered 'single copy' or not
    
    :param family: GeneFamily object
    :param dup_margin: maximal number of genomes in which the gene family can have multiple members and still be considered a 'single copy' gene family

    :return: True if gene family is single copy else False
    """
    nb_multi = 0
    for gene_list in family.get_org_dict().values():
        if len(gene_list) > 1:
            nb_multi += 1
    dup_ratio = nb_multi / family.number_of_organisms
    if dup_ratio < dup_margin:
        return True
    return False


def get_families_to_write(pangenome: Pangenome, partition_filter: str = "core", soft_core: float = 0.95,
                          dup_margin: float = 0.95, single_copy: bool = True) -> Set[GeneFamily]:
    """
    Get families corresponding to the given partition

    :param pangenome: Partitioned pangenome
    :param partition_filter: choice of partition to compute Multiple Sequence Alignement of the gene families
    :param soft_core: Soft core threshold to use
    :param dup_margin: maximal number of genomes in which the gene family can have multiple members and still be considered a 'single copy' gene family
    :param single_copy: Use "single copy" (defined by dup_margin) gene families only

    :return: set of families unique to one partition
    """
    families = set()
    nb_org = pangenome.number_of_organisms

    if partition_filter == "all":
        return pangenome.gene_families
    if partition_filter in ["persistent", "shell", "cloud"]:
        for family in pangenome.gene_families:
            if family.named_partition == partition_filter:
                if single_copy:
                    if is_single_copy(family, dup_margin):
                        families.add(family)
                else:
                    families.add(family)
    elif partition_filter in ["core", "accessory", "softcore"]:
        if partition_filter == "core":
            for family in pangenome.gene_families:
                if family.number_of_organisms == nb_org:
                    if single_copy:
                        if is_single_copy(family, dup_margin):
                            families.add(family)
                    else:
                        families.add(family)
        elif partition_filter == "accessory":
            for family in pangenome.gene_families:
                if family.number_of_organisms < nb_org:
                    if single_copy:
                        if is_single_copy(family, dup_margin):
                            families.add(family)
                    else:
                        families.add(family)
        elif partition_filter == "softcore":
            for family in pangenome.gene_families:
                if family.number_of_organisms >= nb_org * soft_core:
                    if single_copy:
                        if is_single_copy(family, dup_margin):
                            families.add(family)
                    else:
                        families.add(family)
    return families


def translate(seq: str, code: Dict[str, str]) -> str:
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


def write_fasta_families(family: GeneFamily, tmpdir: tempfile.TemporaryDirectory, code_table: Dict[str, str],
                         source: str = 'protein', use_gene_id: bool = False) -> Path:
    """Write fasta files for each gene family

    :param family: gene family to write
    :param tmpdir: path to temporary directory
    :param source: indicates whether to use protein or dna sequences to compute the msa
    :param use_gene_id: Use gene identifiers rather than organism names for sequences in the family MSA
    :param code_table: Genetic code to use

    :return: path to fasta file
    """
    # have a directory for each gene family, to make deletion of tmp files simpler

    f_name = Path(tmpdir.name) / f"{family.name}.fasta"
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


def launch_mafft(fname: Path, output: Path, fam_name: str):
    """
    Compute the MSA with mafft

    :param fname: family gene sequence in fasta
    :param output: directory to save alignment
    :param fam_name: Name of the gene family
    """
    outname = output / f"{fam_name}.aln"
    cmd = ["mafft", "--thread", "1", fname.absolute().as_posix()]
    logging.getLogger("PPanGGOLiN").debug("command: " + " ".join(cmd))
    subprocess.run(cmd, stdout=open(outname, "w"), stderr=subprocess.DEVNULL, check=True)


def launch_multi_mafft(args: List[Tuple[Path, Path, str]]):
    """ Allow to launch mafft in multiprocessing

    :param args: Pack of argument for launch_mafft

    :return: Organism object for pangenome
    """
    launch_mafft(*args)


def compute_msa(families: set, output: Path, tmpdir: Path, cpu: int = 1, source: str = "protein",
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
    logging.getLogger("PPanGGOLiN").info("Preparing input files for MSA...")
    code_table = genetic_codes(str(code))

    for family in tqdm(families, unit="family", disable=disable_bar):
        start_write = time.time()
        fname = write_fasta_families(family, newtmpdir, code_table, source, use_gene_id)
        write_total = write_total + (time.time() - start_write)
        args.append((fname, output, family.name))

    logging.getLogger("PPanGGOLiN").info("Computing the MSA ...")
    bar = tqdm(range(len(families)), unit="family", disable=disable_bar)
    with get_context('fork').Pool(cpu) as p:
        for _ in p.imap_unordered(launch_multi_mafft, args):
            bar.update()
    bar.close()


def write_whole_genome_msa(pangenome: Pangenome, families: set, phylo_name: str, outdir: Path,
                           use_gene_id: bool = False):
    """
    Writes a whole genome msa file for additional phylogenetic analysis

    :param pangenome: Pangenome object
    :param families: Set of families specific to given partition
    :param phylo_name: output file name for phylo alignment
    :param outdir: output directory name for families alignment
    :param use_gene_id: Use gene identifiers rather than organism names for sequences in the family MSA
    """
    phylo_dict = {}
    for org in pangenome.organisms:
        phylo_dict[org.name] = ""
    for fam in families:
        missing_genomes = set(phylo_dict.keys())
        with open(outdir / f"{fam.name}.aln", "r") as fin:
            genome_id = ""
            seq = ""
            curr_len = 0
            curr_phylo_dict = {}

            for line in fin:
                if line.startswith('>'):
                    if genome_id != "":
                        if genome_id not in missing_genomes:
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
        for genome in missing_genomes:
            curr_phylo_dict[genome] = "-" * curr_len

        for key, val in curr_phylo_dict.items():
            phylo_dict[key] += val

    with open(phylo_name, "w") as fout:
        for key, val in phylo_dict.items():
            fout.write(">" + key + "\n")
            fout.write(val + "\n")


def write_msa_files(pangenome: Pangenome, output: Path, cpu: int = 1, partition: str = "core", tmpdir: Path = None,
                    source: str = "protein", soft_core: float = 0.95, phylo: bool = False, use_gene_id: bool = False,
                    translation_table: str = "11", dup_margin: float = 0.95, single_copy: bool = True,
                    force: bool = False, disable_bar: bool = False):
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
    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir

    need_partitions = False
    if partition in ["persistent", "shell", "cloud"]:
        need_partitions = True

    outdir = output / f"msa_{partition}_{source}/"
    mk_outdir(outdir, force=force)

    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_partitions=need_partitions,
                         need_gene_sequences=True, disable_bar=disable_bar)
    logging.getLogger("PPanGGOLiN").info(f"Doing MSA for {partition} families...")
    families = get_families_to_write(pangenome, partition_filter=partition, soft_core=soft_core, dup_margin=dup_margin,
                                     single_copy=single_copy)

    # check that the code is similar than the one used previously, if there is one
    if 'translation_table' in pangenome.parameters["cluster"]:
        if pangenome.parameters["cluster"]["translation_table"] != translation_table:
            logging.getLogger("PPanGGOLiN").warning("The translation table used during clustering "
                                                    f"('{pangenome.parameters['cluster']['translation_table']}') "
                                                    f"is different than the one provided now ('{translation_table}')")
    code = translation_table

    compute_msa(families, outdir, cpu=cpu, tmpdir=tmpdir, source=source, use_gene_id=use_gene_id, code=code,
                disable_bar=disable_bar)
    logging.getLogger("PPanGGOLiN").info(f"Done writing all {partition} MSA in: {outdir}")

    if phylo:
        logging.getLogger("PPanGGOLiN").info("Writing the whole genome msa file")
        if partition == "softcore":
            phylo_name = output / f"{partition}_{soft_core}_genome_alignment.aln"
        else:
            phylo_name = output / f"{partition}_genome_alignment.aln"
        write_whole_genome_msa(pangenome, families, phylo_name, outdir, use_gene_id=use_gene_id)
        logging.getLogger("PPanGGOLiN").info(f"Done writing the {partition} genome alignment in: '{phylo_name}'")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    write_msa_files(pangenome, args.output, cpu=args.cpu, partition=args.partition, tmpdir=args.tmpdir,
                    source=args.source, soft_core=args.soft_core, phylo=args.phylo, use_gene_id=args.use_gene_id,
                    translation_table=args.translation_table, dup_margin=args.dup_margin, single_copy=args.single_copy,
                    force=args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("msa", formatter_class=argparse.RawTextHelpFormatter)
    parser_msa(parser)
    return parser


def parser_msa(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of msa command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="The following arguments are required :")
    required.add_argument('-p', '--pangenome', required=False, type=Path, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=Path,
                          help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title="Optional arguments. Indicating 'all' writes all elements. "
                                               "Writing a partition ('persistent', 'shell', 'cloud', 'core' or "
                                               "'accessory') write the elements associated to said partition.")
    # could make choice to allow customization
    optional.add_argument("--soft_core", required=False, type=restricted_float, default=0.95,
                          help="Soft core threshold to use if 'softcore' partition is chosen")
    optional.add_argument("--dup_margin", required=False, type=restricted_float, default=0.05,
                          help="minimum ratio of genomes in which the family must have multiple genes "
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
                          help="Use gene identifiers rather than genome names for sequences in the family MSA"
                               " (genome names are used by default)")
    optional.add_argument("--translation_table", required=False, default=11, type=int,
                          help="Translation table (genetic code) to use.")
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    optional.add_argument("--tmpdir", required=False, type=str, default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_msa(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())

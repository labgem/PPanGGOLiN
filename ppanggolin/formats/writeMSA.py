#!/usr/bin/env python3

# default libraries
import argparse
import logging
import tempfile
import time
from multiprocessing import get_context
from pathlib import Path
from typing import Dict, Set, List, Tuple

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.genome import Gene
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import (
    mk_outdir,
    restricted_float,
    run_subprocess,
    check_tools_availability,
)
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.genetic_codes import genetic_codes


def get_families_to_write(
    pangenome: Pangenome,
    partition_filter: str = "core",
    soft_core: float = 0.95,
    dup_margin: float = 0.95,
    single_copy: bool = True,
) -> Set[GeneFamily]:
    """
    Get families corresponding to the given partition

    :param pangenome: Partitioned pangenome
    :param partition_filter: choice of partition to compute Multiple Sequence Alignment of the gene families
    :param soft_core: Soft core threshold to use
    :param dup_margin: maximal number of genomes in which the gene family can have multiple members and still be considered a 'single copy' gene family
    :param single_copy: Use "single copy" (defined by dup_margin) gene families only

    :return: set of families unique to one partition
    """
    if partition_filter == "all":
        return set(pangenome.gene_families)
    else:
        families = set()
        nb_org = pangenome.number_of_organisms
        if partition_filter in ["persistent", "shell", "cloud"]:
            for family in pangenome.gene_families:
                if family.named_partition == partition_filter:
                    if single_copy:
                        if family.is_single_copy(
                            dup_margin=dup_margin, exclude_fragment=True
                        ):
                            families.add(family)
                    else:
                        families.add(family)
        elif partition_filter in ["core", "accessory", "softcore"]:
            if partition_filter == "core":
                for family in pangenome.gene_families:
                    if family.number_of_organisms == nb_org:
                        if single_copy:
                            if family.is_single_copy(
                                dup_margin=dup_margin, exclude_fragment=False
                            ):
                                families.add(family)
                        else:
                            families.add(family)
            elif partition_filter == "accessory":
                for family in pangenome.gene_families:
                    if family.number_of_organisms < nb_org:
                        if single_copy:
                            if family.is_single_copy(
                                dup_margin=dup_margin, exclude_fragment=False
                            ):
                                families.add(family)
                        else:
                            families.add(family)
            elif partition_filter == "softcore":
                for family in pangenome.gene_families:
                    if family.number_of_organisms >= nb_org * soft_core:
                        if single_copy:
                            if family.is_single_copy(
                                dup_margin=dup_margin, exclude_fragment=False
                            ):
                                families.add(family)
                        else:
                            families.add(family)
        return families


def translate(gene: Gene, code: Dict[str, Dict[str, str]]) -> Tuple[str, bool]:
    """translates the given dna sequence with the given translation table

    :param gene: given gene
    :param code: translation table corresponding to genetic code to use

    :return: protein sequence
    """
    # code:  https://www.bioinformatics.org/sms/iupac.html
    start_table = code["start_table"]
    table = code["trans_table"]
    mod = len(gene.dna) % 3
    partial = False
    if mod != 0:
        partial = True
        msg = (
            f"Gene {gene.ID} {'' if gene.local_identifier == '' else 'with local identifier ' + gene.local_identifier}"
            f" has a sequence length of {len(gene.dna)} which modulo 3 was different than 0."
        )
        logging.getLogger("PPANGGOLIN").debug(msg)
    protein = start_table[gene.dna[0:3]]
    for i in range(3, len(gene.dna) - mod, 3):
        codon = gene.dna[i : i + 3]
        try:
            protein += table[codon]
        except KeyError:  # codon was not planned for. Probably can't determine it.
            protein += "X"  # X is for unknown
    return protein, partial


def write_fasta_families(
    family: GeneFamily,
    tmpdir: tempfile.TemporaryDirectory,
    code_table: Dict[str, Dict[str, str]],
    source: str = "protein",
    use_gene_id: bool = False,
) -> Tuple[Path, bool]:
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

    # get genes that are present in only one copy for our family in each organism.
    single_copy_genes = []
    for genes in family.get_org_dict().values():
        if len(genes) == 1:
            single_copy_genes.extend(genes)

    with open(f_name, "w") as f_obj:
        partial = False
        for gene in single_copy_genes:
            if use_gene_id:
                f_obj.write(f">{gene.ID}\n")
            else:
                f_obj.write(f">{gene.organism.name}\n")
            if source == "dna":
                f_obj.write(gene.dna + "\n")
            elif source == "protein":
                protein, part = translate(gene, code_table)
                if not partial:
                    partial = part
                f_obj.write(protein + "\n")
            else:
                raise ValueError(
                    f"Unknown sequence source '{source}' provided. Expected 'dna' or 'protein'."
                )

    return f_name, partial


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
    run_subprocess(cmd, outname, msg="mafft failed with the following error:\n")


def launch_multi_mafft(args: List[Tuple[Path, Path, str]]):
    """Allow to launch mafft in multiprocessing

    :param args: Pack of argument for launch_mafft

    :return: Organism object for pangenome
    """
    launch_mafft(*args)


def compute_msa(
    families: Set[GeneFamily],
    output: Path,
    tmpdir: Path,
    cpu: int = 1,
    source: str = "protein",
    use_gene_id: bool = False,
    code: str = "11",
    disable_bar: bool = False,
):
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
    code_table = genetic_codes(code)

    partial = False
    for family in tqdm(families, unit="family", disable=disable_bar):
        start_write = time.time()
        fname, part = write_fasta_families(
            family, newtmpdir, code_table, source, use_gene_id
        )
        if not partial:
            partial = part
        write_total = write_total + (time.time() - start_write)
        args.append((fname, output, family.name))

    if partial:
        logging.getLogger("PPanGGOLiN").warning(
            "Partial gene was found during translation. "
            "Last nucleotides were removed to translate. "
            "Use --verbose 2 to see genes that are partial"
        )
    logging.getLogger("PPanGGOLiN").info("Computing the MSA ...")
    with get_context("fork").Pool(cpu) as p:
        with tqdm(total=len(families), unit="family", disable=disable_bar) as bar:
            for _ in p.imap_unordered(launch_multi_mafft, args):
                bar.update()


def write_whole_genome_msa(
    pangenome: Pangenome,
    families: set,
    phylo_name: Path,
    outdir: Path,
    use_gene_id: bool = False,
):
    """
    Writes a whole genome msa file for additional phylogenetic analysis

    :param pangenome: Pangenome object
    :param families: Set of families specific to given partition
    :param phylo_name: output file name for phylo alignment
    :param outdir: output directory name for families alignment
    :param use_gene_id: Use gene identifiers rather than organism names for sequences in the family MSA
    """

    # sort families by ID, so the gene order is consistent
    families = sorted(families, key=lambda f: f.ID)

    phylo_dict = {}
    for org in pangenome.organisms:
        phylo_dict[org.name] = ""
    for fam in families:
        observed_genomes = set()
        with open(outdir / f"{fam.name}.aln") as fin:
            genome_id = ""
            seq = ""
            curr_len = 0
            curr_phylo_dict = {}

            for line in fin:
                if line.startswith(">"):
                    # Save sequence of previous record
                    if genome_id != "":
                        if genome_id in observed_genomes:
                            # duplicated genes. Replacing them with gaps.
                            curr_phylo_dict[genome_id] = "-" * curr_len
                        else:
                            curr_phylo_dict[genome_id] = seq
                            curr_len = len(seq)
                            observed_genomes.add(genome_id)
                    if use_gene_id:
                        genome_id = pangenome.get_gene(line[1:].strip()).organism.name
                    else:
                        genome_id = line[1:].strip()
                    seq = ""
                else:
                    seq += line.strip()

            # process the final record
            if genome_id != "":
                if genome_id in observed_genomes:
                    # duplicated genes. Replacing them with gaps.
                    curr_phylo_dict[genome_id] = "-" * curr_len
                else:
                    curr_phylo_dict[genome_id] = seq
                    curr_len = len(seq)
                    observed_genomes.add(genome_id)

        # write gaps for all missing genomes
        missing_genomes = [
            g for g in set(phylo_dict.keys()) if g not in observed_genomes
        ]
        for genome in missing_genomes:
            curr_phylo_dict[genome] = "-" * curr_len

        for key, val in curr_phylo_dict.items():
            phylo_dict[key] += val

    with open(phylo_name, "w") as fout:
        for key, val in phylo_dict.items():
            fout.write(">" + key + "\n")
            fout.write(val + "\n")


def write_msa_files(
    pangenome: Pangenome,
    output: Path,
    cpu: int = 1,
    partition: str = "core",
    tmpdir: Path = None,
    source: str = "protein",
    soft_core: float = 0.95,
    phylo: bool = False,
    use_gene_id: bool = False,
    translation_table: str = "11",
    dup_margin: float = 0.95,
    single_copy: bool = True,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Main function to write MSA files

    :param pangenome: Pangenome object with partition
    :param output: Path to output directory
    :param cpu: number of available core
    :param partition: choice of partition to compute Multiple Sequence Alignment of the gene families
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
    check_tools_availability(["mafft"])

    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir

    need_partitions = False
    if partition in ["persistent", "shell", "cloud"]:
        need_partitions = True

    outdir = output / f"msa_{partition}_{source}/"
    mk_outdir(outdir, force=force)

    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_partitions=need_partitions,
        need_gene_sequences=True,
        disable_bar=disable_bar,
    )
    logging.getLogger("PPanGGOLiN").info(f"Doing MSA for {partition} families...")
    families = get_families_to_write(
        pangenome,
        partition_filter=partition,
        soft_core=soft_core,
        dup_margin=dup_margin,
        single_copy=single_copy,
    )

    # check that the code is similar than the one used previously, if there is one
    if "translation_table" in pangenome.parameters["cluster"]:
        if pangenome.parameters["cluster"]["translation_table"] != str(
            translation_table
        ):
            logging.getLogger("PPanGGOLiN").warning(
                "The translation table used during clustering "
                f"('{pangenome.parameters['cluster']['translation_table']}') "
                f"is different than the one provided now ('{translation_table}')"
            )

    compute_msa(
        families,
        outdir,
        cpu=cpu,
        tmpdir=tmpdir,
        source=source,
        use_gene_id=use_gene_id,
        code=str(translation_table),
        disable_bar=disable_bar,
    )
    logging.getLogger("PPanGGOLiN").info(
        f"Done writing all {partition} MSA in: {outdir}"
    )

    if phylo:
        logging.getLogger("PPanGGOLiN").info("Writing the whole genome msa file")
        if partition == "softcore":
            phylo_name = output / f"{partition}_{soft_core}_genome_alignment.aln"
        else:
            phylo_name = output / f"{partition}_genome_alignment.aln"
        write_whole_genome_msa(
            pangenome, families, phylo_name, outdir, use_gene_id=use_gene_id
        )
        logging.getLogger("PPanGGOLiN").info(
            f"Done writing the {partition} genome alignment in: '{phylo_name}'"
        )


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    mk_outdir(args.output, args.force)
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    write_msa_files(
        pangenome,
        args.output,
        cpu=args.cpu,
        partition=args.partition,
        tmpdir=args.tmpdir,
        source=args.source,
        soft_core=args.soft_core,
        phylo=args.phylo,
        use_gene_id=args.use_gene_id,
        translation_table=args.translation_table,
        dup_margin=args.dup_margin,
        single_copy=args.single_copy,
        force=args.force,
        disable_bar=args.disable_prog_bar,
    )


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
    required = parser.add_argument_group(
        title="Required arguments", description="The following arguments are required :"
    )
    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome .h5 file"
    )
    required.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output directory where the file(s) will be written",
    )

    optional = parser.add_argument_group(
        title="Optional arguments. Indicating 'all' writes all elements. "
        "Writing a partition ('persistent', 'shell', 'cloud', 'core' or "
        "'accessory') write the elements associated to said partition."
    )
    # could make choice to allow customization
    optional.add_argument(
        "--soft_core",
        required=False,
        type=restricted_float,
        default=0.95,
        help="Soft core threshold to use if 'softcore' partition is chosen",
    )
    optional.add_argument(
        "--dup_margin",
        required=False,
        type=restricted_float,
        default=0.05,
        help="minimum ratio of genomes in which the family must have multiple genes "
        "for it to be considered 'duplicated'",
    )
    optional.add_argument(
        "--single_copy",
        required=False,
        action="store_true",
        default=False,
        help="Use report gene families that are considered 'single copy', for details see "
        "option --dup_margin",
    )
    optional.add_argument(
        "--partition",
        required=False,
        default="core",
        choices=[
            "all",
            "persistent",
            "shell",
            "cloud",
            "core",
            "accessory",
            "softcore",
        ],
        help="compute Multiple Sequence Alignment of the gene families in the given partition",
    )
    optional.add_argument(
        "--source",
        required=False,
        default="protein",
        choices=["dna", "protein"],
        help="indicates whether to use protein or dna sequences to compute the msa",
    )
    optional.add_argument(
        "--phylo",
        required=False,
        action="store_true",
        help="Writes a whole genome msa file for additional phylogenetic analysis",
    )
    optional.add_argument(
        "--use_gene_id",
        required=False,
        action="store_true",
        help="Use gene identifiers rather than genome names for sequences in the family MSA"
        " (genome names are used by default)",
    )
    optional.add_argument(
        "--translation_table",
        required=False,
        default=11,
        type=int,
        help="Translation table (genetic code) to use.",
    )
    optional.add_argument(
        "-c",
        "--cpu",
        required=False,
        default=1,
        type=int,
        help="Number of available cpus",
    )
    optional.add_argument(
        "--tmpdir",
        required=False,
        type=str,
        default=Path(tempfile.gettempdir()),
        help="directory for storing temporary files",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_msa(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())

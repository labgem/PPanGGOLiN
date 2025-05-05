#!/usr/bin/env python3

# default libraries
import logging
import argparse
from pathlib import Path
from typing import Set, Iterable

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.genome import Organism, Contig
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region
from ppanggolin.formats import check_pangenome_info, write_pangenome, erase_pangenome
from ppanggolin.utils import restricted_float


class MatriceNode:
    def __init__(self, state, score, prev, gene):
        self.state = state  # state of the node. 1 for RGP and 0 for not RGP.
        self.score = score if score > 0 else 0  # current score of the node
        self.prev = prev  # previous matriceNode
        self.gene = gene  # gene this node corresponds to

    def changes(self, score):
        # state of the node. 1 for RGP and 0 for not RGP.
        self.state = 1 if score >= 0 else 0
        # current score of the node. If the given score is negative, set to 0.
        self.score = score if score >= 0 else 0


def extract_rgp(contig, node, rgp_id, naming) -> Region:
    """
    Extract the region from the given starting node
    """
    new_region = None
    if naming == "contig":
        new_region = Region(contig.name + "_RGP_" + str(rgp_id))
    elif naming == "organism":
        new_region = Region(
            node.gene.organism.name + "_" + contig.name + "_RGP_" + str(rgp_id)
        )
    while node.state:
        new_region.add(node.gene)
        node.state = 0
        node.score = 0
        node = node.prev
        if node is None:  # it's the end of the contig and the end of the region.
            break
    return new_region


def rewrite_matrix(contig, matrix, index, persistent, continuity, multi):
    """
    ReWrite the matrice from the given index of the node that started a region.
    """
    prev = matrix[index]
    index += 1
    if index > len(matrix) and contig.is_circular:
        index = 0
    # else the node was the last one of the contig, and there is nothing to do
    if index < len(matrix):
        next_node = matrix[index]
        nb_perc = 0
        while next_node.state:  # while the old state is not 0, recompute the scores.
            if (
                next_node.gene.family.named_partition == "persistent"
                and next_node.gene.family not in multi
            ):
                modif = -pow(persistent, nb_perc)
                nb_perc += 1
            else:
                modif = continuity
                nb_perc = 0

            curr_score = modif + prev.score
            # scores can't be negative. If they are, they'll be set to 0.
            matrix[index].changes(curr_score)
            index += 1
            if index >= len(matrix):
                if contig.is_circular:
                    index = 0
                else:
                    # else we're at the end of the contig, so there are no more computations. Get out of the loop
                    break

            prev = next_node
            next_node = matrix[index]


def init_matrices(
    contig: Contig, multi: set, persistent_penalty: int = 3, variable_gain: int = 1
) -> list:
    """
    Initialize the vector of score/state nodes

    :param contig: Current contig from one organism
    :param persistent_penalty: Penalty score to apply to persistent genes
    :param variable_gain: Gain score to apply to variable genes
    :param multi: multigenic persistent families of the pangenome graph.

    :return: Initialized matrice
    """
    mat = []
    prev = None
    nb_perc = 0
    zero_ind = None
    curr_state = None
    for gene in contig.genes:
        if gene.family.named_partition == "persistent" and gene.family not in multi:
            modif = -pow(persistent_penalty, nb_perc)
            nb_perc += 1
        else:
            modif = variable_gain
            nb_perc = 0

        curr_score = modif + prev.score if prev is not None else modif
        if curr_score >= 0:
            curr_state = 1
        else:
            curr_state = 0
            zero_ind = True
        prev = MatriceNode(curr_state, curr_score, prev, gene)
        if prev.state == 0:
            zero_ind = prev
        mat.append(prev)
        logging.getLogger("PPanGGOLiN").debug(
            f"gene:{gene.ID};zero_ind:{zero_ind};curr_state:{curr_state};curr_score:{curr_score}."
        )

    if zero_ind is None:
        zero_ind = prev  # don't go further than the current node, if no node were at 0.

    # if the contig is circular, and we're in a rgp state,
    # we need to continue from the "starting" gene until we leave rgp state.
    if contig.is_circular and curr_state and zero_ind is not None:
        # the previous node of the first processed gene is the last node.
        mat[0].prev = prev
        c = 0
        nb_perc = 0
        while curr_state:  # while state is rgp.
            mat_node = mat[c]
            if mat_node == zero_ind:
                # then we've parsed the entire contig twice.
                logging.getLogger("PPanGGOLiN").debug(
                    f"{contig.name} was parsed entirely twice."
                )
                # The whole sequence is a rgp, so we're stopping the iteration now, otherwise we'll loop indefinitely
                break

            if (
                mat_node.gene.family.named_partition == "persistent"
                and mat_node.gene.family not in multi
            ):
                modif = -pow(persistent_penalty, nb_perc)
                nb_perc += 1
            else:
                modif = variable_gain
                nb_perc = 0

            curr_score = modif + mat_node.prev.score
            curr_state = 1 if curr_score >= 0 else 0
            mat_node.changes(curr_score)
            logging.getLogger("PPanGGOLiN").debug(
                f"gene:{mat_node.gene.ID};curr_state:{curr_state};curr_score:{curr_score}."
            )
            c += 1
    return mat


def mk_regions(
    contig: Contig,
    matrix: list,
    multi: set,
    min_length: int = 3000,
    min_score: int = 4,
    persistent: int = 3,
    continuity: int = 1,
    naming: str = "contig",
) -> Set[Region]:
    """
    Processing matrix and 'emptying' it to get the regions.

    :param contig: Current contig from one organism
    :param matrix: Initialized matrix
    :param multi: multigenic persistent families of the pangenome graph.
    :param min_length: Minimum length (bp) of a region to be considered RGP
    :param min_score: Minimal score wanted for considering a region as being RGP
    :param persistent: Penalty score to apply to persistent genes
    :param continuity: Gain score to apply to variable genes
    :param naming:

    :return:
    """

    def max_index_node(lst):
        """gets the last node with the highest score from a list of matriceNode"""
        if isinstance(lst, list):
            # init with the first element of the list
            max_score = lst[0].score
            max_index = 0
            for idx, node in enumerate(lst):
                if node.score >= max_score:
                    max_score = node.score
                    max_index = idx
            return max_score, max_index
        else:
            raise TypeError(
                f"List of matriceNode is expected. The detected type was {type(lst)}"
            )

    contig_regions = set()
    val, index = max_index_node(matrix)
    while val >= min_score:
        new_region = extract_rgp(contig, matrix[index], len(contig_regions), naming)
        new_region.score = val
        if new_region.length > min_length:
            contig_regions.add(new_region)
        else:
            # Remove region reference in genes to not consider them when iterating genome RGP
            for gene in new_region.genes:
                gene._RGP = None

        rewrite_matrix(contig, matrix, index, persistent, continuity, multi)
        val, index = max_index_node(matrix)
    return contig_regions


def compute_org_rgp(
    organism: Organism,
    multigenics: set,
    persistent_penalty: int = 3,
    variable_gain: int = 1,
    min_length: int = 3000,
    min_score: int = 4,
    naming: str = "contig",
    disable_bar: bool = True,
) -> set:
    """
    Compute regions of genomic plasticity (RGP) on the given organism based on the provided parameters.

    :param organism: The Organism object representing the organism.
    :param multigenics: A set of multigenic persistent families of the pangenome graph.
    :param persistent_penalty: Penalty score to apply to persistent multigenic families (default: 3).
    :param variable_gain: Gain score to apply to variable multigenic families (default: 1).
    :param min_length: Minimum length threshold (in base pairs) for the regions to be considered RGP (default: 3000).
    :param min_score: Minimum score threshold for considering a region as RGP (default: 4).
    :param naming: Naming scheme for the regions, either "contig" or "organism" (default: "contig").
    :param disable_bar: Whether to disable the progress bar. It is recommended to disable it when calling this function in a loop on multiple organisms (default: True).

    :return: A set of RGPs of the provided organism.
    """
    org_regions = set()
    for contig in tqdm(
        organism.contigs,
        total=organism.number_of_contigs,
        unit="contig",
        disable=disable_bar,
    ):
        if contig.number_of_genes != 0:  # some contigs have no coding genes...
            # can definitely multiprocess this part, as not THAT much information is needed...
            matrix = init_matrices(
                contig, multigenics, persistent_penalty, variable_gain
            )
            org_regions |= mk_regions(
                contig,
                matrix,
                multigenics,
                min_length,
                min_score,
                persistent_penalty,
                variable_gain,
                naming=naming,
            )
    return org_regions


def naming_scheme(organisms: Iterable[Organism]) -> str:
    """
    Determine the naming scheme for the contigs in the pangenome.

    :param organisms: Iterable of organims objects
    :return: Naming scheme for the contigs ("contig" or "organism").
    """
    contigsids = set()
    for org in organisms:
        for contig in org.contigs:
            oldlen = len(contigsids)
            contigsids.add(contig.name)
            if oldlen == len(contigsids):
                logging.getLogger("PPanGGOLiN").warning(
                    "You have contigs with identical identifiers in your "
                    "assemblies. Identifiers will be supplemented with your "
                    "provided organism names."
                )
                return "organism"
    return "contig"


def check_pangenome_former_rgp(pangenome: Pangenome, force: bool = False):
    """checks pangenome status and .h5 files for former rgp, delete them if allowed or raise an error

    :param pangenome: Pangenome object
    :param force: Allow to force write on Pangenome file
    """
    if pangenome.status["predictedRGP"] == "inFile" and not force:
        raise Exception(
            "You are trying to predict RGPs in a pangenome that already have them predicted. "
            "If you REALLY want to do that, use --force "
            "(it will erase RGPs and every feature computed from them)."
        )
    elif pangenome.status["predictedRGP"] == "inFile" and force:
        erase_pangenome(pangenome, rgp=True)


def predict_rgp(
    pangenome: Pangenome,
    persistent_penalty: int = 3,
    variable_gain: int = 1,
    min_length: int = 3000,
    min_score: int = 4,
    dup_margin: float = 0.05,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Main function to predict region of genomic plasticity

    :param pangenome: blank pangenome object
    :param persistent_penalty: Penalty score to apply to persistent genes
    :param variable_gain: Gain score to apply to variable genes
    :param min_length: Minimum length (bp) of a region to be considered RGP
    :param min_score: Minimal score wanted for considering a region as being RGP
    :param dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    # check statuses and load info
    check_pangenome_former_rgp(pangenome, force)
    check_pangenome_info(
        pangenome,
        need_annotations=True,
        need_families=True,
        need_graph=False,
        need_partitions=True,
        disable_bar=disable_bar,
    )

    logging.getLogger("PPanGGOLiN").info("Detecting multigenic families...")
    multigenics = pangenome.get_multigenics(dup_margin)
    logging.getLogger("PPanGGOLiN").info("Compute Regions of Genomic Plasticity ...")
    name_scheme = naming_scheme(pangenome.organisms)
    for org in tqdm(
        pangenome.organisms,
        total=pangenome.number_of_organisms,
        unit="genomes",
        disable=disable_bar,
    ):
        for region in compute_org_rgp(
            org,
            multigenics,
            persistent_penalty,
            variable_gain,
            min_length,
            min_score,
            naming=name_scheme,
        ):
            pangenome.add_region(region)
    logging.getLogger("PPanGGOLiN").info(f"Predicted {pangenome.number_of_rgp} RGP")

    # save parameters and save status
    pangenome.parameters["rgp"] = {}
    pangenome.parameters["rgp"]["persistent_penalty"] = persistent_penalty
    pangenome.parameters["rgp"]["variable_gain"] = variable_gain
    pangenome.parameters["rgp"]["min_length"] = min_length
    pangenome.parameters["rgp"]["min_score"] = min_score
    pangenome.parameters["rgp"]["dup_margin"] = dup_margin
    pangenome.status["predictedRGP"] = "Computed"


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    predict_rgp(
        pangenome,
        persistent_penalty=args.persistent_penalty,
        variable_gain=args.variable_gain,
        min_length=args.min_length,
        min_score=args.min_score,
        dup_margin=args.dup_margin,
        force=args.force,
        disable_bar=args.disable_prog_bar,
    )
    write_pangenome(
        pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("rgp", formatter_class=argparse.RawTextHelpFormatter)
    parser_rgp(parser)
    return parser


def parser_rgp(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of rgp command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="One of the following arguments is required :",
    )
    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome .h5 file"
    )

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "--persistent_penalty",
        required=False,
        type=int,
        default=3,
        help="Penalty score to apply to persistent genes",
    )
    optional.add_argument(
        "--variable_gain",
        required=False,
        type=int,
        default=1,
        help="Gain score to apply to variable genes",
    )
    optional.add_argument(
        "--min_score",
        required=False,
        type=int,
        default=4,
        help="Minimal score wanted for considering a region as being a RGP",
    )
    optional.add_argument(
        "--min_length",
        required=False,
        type=int,
        default=3000,
        help="Minimum length (bp) of a region to be considered a RGP",
    )
    optional.add_argument(
        "--dup_margin",
        required=False,
        type=restricted_float,
        default=0.05,
        help="Minimum ratio of genomes where the family is present in which the family must "
        "have multiple genes for it to be considered 'duplicated'",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_rgp(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())

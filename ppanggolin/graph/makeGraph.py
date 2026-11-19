#!/usr/bin/env python3

# default libraries
import logging
import argparse
from pathlib import Path

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.genome import Contig, Gene
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import read_pangenome, write_pangenome, erase_pangenome


def check_pangenome_former_graph(pangenome: Pangenome, force: bool = False):
    """
    Checks pangenome status and .h5 files for former neighbors graph, delete it if allowed or raise an error

    :param pangenome: Pangenome object
    :param force: Allow to force write on Pangenome file
    """
    if pangenome.status["neighborsGraph"] == "inFile" and not force:
        raise AttributeError(
            "You are trying to make a neighbors graph that is already built. "
            "If you REALLY want to do that, use --force "
            "(it will erase everything except annotation data !)"
        )
    elif pangenome.status["neighborsGraph"] == "inFile" and force:
        erase_pangenome(pangenome, graph=True)


def check_pangenome_for_neighbors_graph(pangenome, force, disable_bar=False):
    """
    Checks and read the pangenome for neighbors graph computing.

    :param pangenome: Pangenome object
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    check_pangenome_former_graph(pangenome, force)
    # TODO Check if possible to change for check_pangenome_info
    if pangenome.status["genomesAnnotated"] in [
        "Computed",
        "Loaded",
    ] and pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        pass  # nothing to do, can just continue.
    elif (
        pangenome.status["genomesAnnotated"] == "inFile"
        and pangenome.status["genesClustered"] == "inFile"
    ):
        read_pangenome(
            pangenome, annotation=True, gene_families=True, disable_bar=disable_bar
        )
    elif pangenome.status["genesClustered"] == "No" and pangenome.status[
        "genomesAnnotated"
    ] in ["inFile", "Computed", "Loaded"]:
        raise Exception(
            "You did not cluster the genes. See the 'ppanggolin cluster' if you want to do that."
        )
    else:
        # You probably can use readPangenome anyway.
        msg = (
            "Dev : You are probably writing a new workflow with a combination that I did not test."
            " You can probably use readPangenome instead of raising this Error. "
            "However please test it carefully.\n"
        )
        msg += (
            " User : I have no idea how you got there. You probably did something unexpected. "
            "Post an issue with what you did at https://github.com/labgem/PPanGGOLiN\n"
        )
        raise NotImplementedError(msg)


def remove_high_copy_number(pangenome, number):
    """Removes families present more than 'number' times from the pangenome graph

    :param pangenome: Pangenome object
    :param number: Maximum authorized repeat presence
    """
    for fam in pangenome.gene_families:
        for gene_list in fam.get_org_dict().values():
            if len(gene_list) >= number:
                fam.removed = True


def _is_fragment_self_link(gene1: Gene, gene2: Gene) -> bool:
    """Tell whether two genes are pieces of the same fragmented gene.

    Two neighboring genes of the same family, one of them flagged as a fragment, are not
    linked in the graph: they likely come from a gene split in pieces, not from two neighbor genes.

    :param gene1: First gene
    :param gene2: Second gene
    :return: True if the two genes must not be linked
    """
    return gene1.family == gene2.family and (gene1.is_fragment or gene2.is_fragment)


def add_contig_edges(pangenome: Pangenome, contig: Contig):
    """Add to the neighbors graph the edges between the consecutive genes of a contig.

    Genes of removed families are skipped, so their neighbors are linked together, and
    two genes of a fragmented gene are not linked. On a circular contig, the last and
    the first kept genes are linked under the same rules.

    :param pangenome: Pangenome receiving the edges
    :param contig: Contig whose genes are linked
    """
    first = prev = None
    for gene in contig.genes:
        if gene.family is None:
            raise AttributeError("a Gene does not have a GeneFamily object associated")
        if gene.family.removed:
            continue
        if prev is None:
            first = gene
        elif not _is_fragment_self_link(prev, gene):
            pangenome.add_edge(gene, prev)
        prev = gene

    # first is None when the contig only has genes of removed families,
    # and first is prev when a single gene is kept: there is nothing to close then.
    if (
        contig.is_circular
        and first is not None
        and first is not prev
        and not _is_fragment_self_link(prev, first)
    ):
        pangenome.add_edge(first, prev)


def compute_neighbors_graph(
    pangenome: Pangenome,
    remove_copy_number: int = 0,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Creates the Pangenome Graph. Will either load the information from the pangenome file if they are not loaded,
    or use the information loaded if they are.

    :param pangenome: Pangenome object
    :param remove_copy_number: Maximum authorized repeat presence of gene families. if zero no remove
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    check_pangenome_for_neighbors_graph(pangenome, force, disable_bar=disable_bar)

    if remove_copy_number > 0:
        remove_high_copy_number(pangenome, remove_copy_number)

    logging.getLogger("PPanGGOLiN").info("Computing the neighbors graph...")
    bar = tqdm(
        pangenome.organisms,
        total=pangenome.number_of_organisms,
        unit="genome",
        disable=disable_bar,
    )
    for org in bar:
        bar.set_description(f"Processing {org.name}")
        bar.refresh()
        for contig in org.contigs:
            add_contig_edges(pangenome, contig)
    logging.getLogger("PPanGGOLiN").info("Done making the neighbors graph.")
    pangenome.status["neighborsGraph"] = "Computed"

    pangenome.parameters["graph"] = {}
    if remove_copy_number > 0:
        pangenome.parameters["graph"]["remove_high_copy_number"] = remove_copy_number


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    compute_neighbors_graph(
        pangenome,
        args.remove_high_copy_number,
        args.force,
        disable_bar=args.disable_prog_bar,
    )
    write_pangenome(
        pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for graph command

    :return : parser arguments for graph command
    """
    parser = sub_parser.add_parser(
        "graph", formatter_class=argparse.RawTextHelpFormatter
    )
    parser.description = "Create the pangenome graph"
    parser.category = "Expert"
    parser_graph(parser)
    return parser


def parser_graph(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of graph command

    :param parser: parser for graph argument
    """
    required = parser.add_argument_group(
        title="Required arguments", description="Following arguments is required:"
    )
    required.add_argument(
        "-p",
        "--pangenome",
        required=False,
        type=Path,
        help="Path to the pangenome .h5 file.",
    )
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "-r",
        "--remove_high_copy_number",
        type=int,
        default=0,
        help="Positive Number: Remove families having a number of copy of gene in a single genome "
        "above or equal to this threshold in at least one genome "
        "(0 or negative values are ignored).",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_graph(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())

#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import argparse

# installed libraries
from tqdm import tqdm

# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import read_pangenome, write_pangenome, erase_pangenome


def check_pangenome_former_graph(pangenome: Pangenome, force: bool = False):
    """
    Checks pangenome status and .h5 files for former neighbors graph, delete it if allowed or raise an error

    :param pangenome: Pangenome object
    :param force: Allow to force write on Pangenome file
    """
    if pangenome.status["neighborsGraph"] == "inFile" and not force:
        raise Exception("You are trying to make a neighbors graph that is already built. "
                        "If you REALLY want to do that, use --force (it will erase everything except annotation data !)"
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
    if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] and \
            pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        pass  # nothing to do, can just continue.
    elif pangenome.status["genomesAnnotated"] == "inFile" and pangenome.status["genesClustered"] == "inFile":
        read_pangenome(pangenome, annotation=True, gene_families=True, disable_bar=disable_bar)
    elif pangenome.status["genesClustered"] == "No" and \
            pangenome.status["genomesAnnotated"] in ['inFile', 'Computed', 'Loaded']:
        raise Exception("You did not cluster the genes. See the 'ppanggolin cluster' if you want to do that.")
    else:
        # You probably can use readPangenome anyway.
        msg = "Dev : You are probably writing a new workflow with a combination that I did not test." \
              " You can probably use readPangenome instead of raising this Error. " \
              "However please test it carefully.\n"
        msg += " User : I have no idea how you got there. You probably did something unexpected. " \
               "Post an issue with what you did at https://github.com/labgem/PPanGGOLiN\n"
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


def compute_neighbors_graph(pangenome: Pangenome, remove_copy_number: int = 0,
                            force: bool = False, disable_bar: bool = False):
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

    logging.getLogger().info("Computing the neighbors graph...")
    bar = tqdm(pangenome.organisms, total=len(pangenome.organisms), unit="organism", disable=disable_bar)
    for org in bar:
        bar.set_description(f"Processing {org.name}")
        bar.refresh()
        for contig in org.contigs:
            prev = None
            for gene in contig.genes:
                try:
                    if not gene.family.removed:
                        if prev is not None and not (prev.family == gene.family and (prev.is_fragment or
                                                                                     gene.is_fragment)):
                            pangenome.add_edge(gene, prev)
                        prev = gene
                except AttributeError:
                    raise AttributeError("a Gene does not have a GeneFamily object associated")
                except Exception:
                    raise Exception("Unexpected error. Please report on our github.")
            if prev is not None and contig.is_circular and len(contig.genes) > 0:
                # if prev is None, the contig is entirely made of duplicated genes, so no edges are added
                pangenome.add_edge(contig.genes[0], prev)
    logging.getLogger().info("Done making the neighbors graph.")
    pangenome.status["neighborsGraph"] = "Computed"

    pangenome.parameters["graph"] = {}
    pangenome.parameters["graph"]["removed_high_copy_number_families"] = False
    if remove_copy_number > 0:
        pangenome.parameters["graph"]["removed_high_copy_number_families"] = True
        pangenome.parameters["graph"]["removed_high_copy_number_of_families_above"] = remove_copy_number


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    compute_neighbors_graph(pangenome, args.remove_high_copy_number, args.force, disable_bar=args.disable_prog_bar)
    write_pangenome(pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for graph command

    :return : parser arguments for graph command
    """
    parser = sub_parser.add_parser("graph", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_graph(parser)
    return parser


def parser_graph(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of graph command

    :param parser: parser for graph argument
    """
    parser.add_argument('-p', '--pangenome', required=False, type=str, help="The pangenome .h5 file")
    parser.add_argument('-r', '--remove_high_copy_number', type=int, default=0,
                        help="Positive Number: Remove families having a number of copy of gene in a single organism "
                             "above or equal to this threshold in at least one organism "
                             "(0 or negative values are ignored).")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    from ppanggolin.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_graph(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('-f', '--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())

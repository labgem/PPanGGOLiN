#!/usr/bin/env python3

# default libraries
import argparse
from pathlib import Path

# installed libraries
import tables
import yaml

# local libraries
from ppanggolin.formats import read_info, read_parameters


def print_yaml(yaml_dict: dict) -> None:
    yaml_output = yaml.dump(
        yaml_dict, default_flow_style=False, sort_keys=False, indent=4
    )
    print(yaml_output)


def read_status(h5f: tables.File):
    """
    Read status on what have been computed in the pangenome file

    :param h5f: the h5f file object of the pangenome file
    """
    status_group = h5f.root.status

    status_to_print = {
        "Genomes_Annotated": status_group._v_attrs.genomesAnnotated,
        "Genes_Clustered": status_group._v_attrs.genesClustered,
        "Genes_with_Sequences": status_group._v_attrs.geneSequences,
        "Gene_Families_with_Sequences": status_group._v_attrs.geneFamilySequences,
        "Neighbors_Graph": status_group._v_attrs.NeighborsGraph,
        "Pangenome_Partitioned": status_group._v_attrs.Partitioned,
        "RGP_Predicted": status_group._v_attrs.predictedRGP,
        "Spots_Predicted": status_group._v_attrs.predictedRGP,
        # Please confirm if this should be different from "RGP Predicted"
        "Modules_Predicted": status_group._v_attrs.modules,
    }
    status_to_print = {key: bool(val) for key, val in status_to_print.items()}

    status_to_print["PPanGGOLiN_Version"] = str(status_group._v_attrs.version)

    return {"Status": status_to_print}


def read_metadata_status(h5f: tables.File):
    """
    Print which object has metadata and the source of the metadata

    :param h5f: the h5f file object of the pangenome file
    """
    status_group = h5f.root.status
    metadata_info = None

    if hasattr(status_group._v_attrs, "metadata") and status_group._v_attrs.metadata:
        metastatus = status_group.metastatus
        metasources = status_group.metasources
        metadata_info = {
            attr: ", ".join(metasources._v_attrs[attr])
            for attr in metastatus._v_attrs._f_list()
        }

    return {"Metadata": metadata_info}


def print_info(
    pangenome: str,
    status: bool = False,
    content: bool = False,
    parameters: bool = False,
    metadata: bool = False,
):
    """
    Main function to return information about pangenome

    :param pangenome: Pangenome file
    :param status: Get pangenome status
    :param content: Get pangenome content
    :param parameters: Get pangenome parameters
    """
    if not (status or content or parameters or metadata):
        status, content, parameters, metadata = (True, True, True, True)

    h5f = tables.open_file(pangenome, "r+")
    if status:
        print_yaml(read_status(h5f))
    if content:
        print_yaml(read_info(h5f))
    if parameters:
        read_parameters(h5f)
    if metadata:
        print_yaml(read_metadata_status(h5f))
    h5f.close()


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    print_info(
        args.pangenome, args.status, args.content, args.parameters, args.metadata
    )


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for info command

    :return : parser arguments for info command
    """
    parser = sub_parser.add_parser(
        "info", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_info(parser)
    return parser


def parser_info(parser: argparse.ArgumentParser):
    """
    Parser for the specific argument of the 'info' command.

    :param parser: Parser for the 'info' argument.
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="Specify the following required argument:",
    )
    required.add_argument(
        "-p",
        "--pangenome",
        required=True,
        type=Path,
        help="Path to the pangenome .h5 file",
    )

    options = parser.add_argument_group(
        title="Information Display Options (default: all)"
    )
    options.add_argument(
        "-a",
        "--parameters",
        required=False,
        action="store_true",
        help="Display the parameters used or computed for each step of pangenome generation",
    )
    options.add_argument(
        "-c",
        "--content",
        required=False,
        action="store_true",
        help="Display detailed information about the pangenome's content",
    )
    options.add_argument(
        "-s",
        "--status",
        required=False,
        action="store_true",
        help="Display information about the statuses of different elements in the pangenome, "
        "indicating what has been computed or not",
    )
    options.add_argument(
        "-m",
        "--metadata",
        required=False,
        action="store_true",
        help="Display a summary of the metadata saved in the pangenome",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_info(main_parser)
    launch(main_parser.parse_args())

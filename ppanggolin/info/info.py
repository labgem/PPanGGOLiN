#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse

# installed libraries
import tables

# local libraries
from ppanggolin.formats import read_info, read_parameters


def print_info(pangenome, status=False, content=False, parameters=False):
    if status or content or parameters:
        h5f = tables.open_file(pangenome, "r")
        if status:
            status_group = h5f.root.status
            print(f"genomes annotated : {'true' if status_group._v_attrs.genomesAnnotated else 'false'}")
            print(f"genes clustered : {'true' if status_group._v_attrs.genesClustered else 'false'}")
            print(f"genes have their sequences : {'true' if status_group._v_attrs.geneSequences else 'false'}")
            print(f"gene families have their sequences : "
                  f"{'true' if status_group._v_attrs.geneFamilySequences else 'false'}")
            print(f"neighbors graph : {'true' if status_group._v_attrs.NeighborsGraph else 'false'}")
            if 'Partitionned' in status_group._v_attrs._f_list():
                # Partitionned keep working with older version
                if status_group._v_attrs.Partitionned:
                    status_group._v_attrs.Partitioned = True
                del status_group._v_attrs.Partitionned
            if status_group._v_attrs.Partitioned:
                print("pangenome partitioned : true")
            else:
                print("pangenome partitioned : false")
            if hasattr(status_group._v_attrs, "predictedRGP"):
                print(f"RGP predicted : {'true' if status_group._v_attrs.predictedRGP else 'false'}")

            if hasattr(status_group._v_attrs, "spots"):
                print(f"Spots predicted : {'true' if status_group._v_attrs.spots else 'false'}")

            if hasattr(status_group._v_attrs, "modules"):
                print(f"Modules predicted : {'true' if status_group._v_attrs.modules else 'false'}")

            if hasattr(status_group._v_attrs, "version"):
                print(f"PPanGGOLiN version : {status_group._v_attrs.version}")

        if content:
            read_info(h5f)
        if parameters:
            read_parameters(h5f)
        h5f.close()
    else:
        print("Please select what information you want by using --parameters, --content or --status")


def launch(args):
    print_info(args.pangenome, args.status, args.content, args.parameters)


def subparser(sub_parser):
    parser = sub_parser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_info(parser)
    return parser


def parser_info(parser):
    required = parser.add_argument_group(title="Required arguments",
                                         description="The following arguments is required :")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")

    options = parser.add_argument_group(title="optional arguments")
    options.add_argument("--parameters", required=False, action="store_true",
                         help="Shows the parameters used (or computed) for each step of the pangenome generation")
    options.add_argument("--content", required=False, action="store_true",
                         help="Shows detailled informations about the pan's content")
    options.add_argument("--status", required=False, action="store_true",
                         help="Shows informations about the statuses of the different elements of the pangenome "
                              "(what has been computed, or not)")


if __name__ == '__main__':
    """To test local change and allow using debugger"""
    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_info(main_parser)
    launch(main_parser.parse_args())

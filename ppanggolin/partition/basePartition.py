#!/usr/bin/env python3
#coding:utf-8

import argparse

def workflow_partition():
    """
        launch partition from a preceding step in ppanggolin
    """
    pass

def launch(arguments):
    """
        main code when launch partition from the command line.
    """
    pass


def partitionSubparser(subparser):
    parser = subparser.add_parser("partition",help = "Partition the pangenome graph")

    return parser
#!/usr/bin/env python3
#coding:utf-8

#default libraries
import sys
if sys.version_info < (3, 6):#minimum is python3.6
    raise AssertionError("Minimum python version to run PPanGGOLiN is 3.6. Your current python version is " + ".".join(map(str,sys.version_info)))
import argparse
import logging
import resource
import pkg_resources
import tempfile
import os

#local modules
import ppanggolin.pangenome
import ppanggolin.nem.partition
import ppanggolin.nem.rarefaction
import ppanggolin.graph
import ppanggolin.annotate
import ppanggolin.cluster
import ppanggolin.workflow.workflow
import ppanggolin.workflow.panRGP
import ppanggolin.figures
import ppanggolin.formats
import ppanggolin.info
import ppanggolin.align
import ppanggolin.RGP

def checkTsvSanity(tsv):
    f = open(tsv,"r")
    nameSet = set()
    duplicatedNames = set()
    nonExistingFiles = set()
    for line in f:
        elements = [el.strip() for el in line.split("\t")]
        if len(elements)<=1:
            raise Exception(f"No tabulation separator found in given file: {tsv}")
        if " " in elements[0]:
            raise Exception(f"Your genome names contain spaces (The first encountered genome name that had this string : '{elements[0]}'). To ensure compatibility with all of the dependencies of PPanGGOLiN this is not allowed. Please remove spaces from your genome names.")
        oldLen = len(nameSet)
        nameSet.add(elements[0])
        if len(nameSet) == oldLen:
            duplicatedNames.add(elements[0])
        if not os.path.exists(elements[1]):
            nonExistingFiles.add(elements[1])
    if len(nonExistingFiles) != 0:
        raise Exception(f"Some of the given files do not exist. The non-existing files are the following : '{' '.join(nonExistingFiles)}'")
    if len(duplicatedNames) != 0:
        raise Exception(f"Some of your genomes have identical names. The duplicated names are the following : '{' '.join(duplicatedNames)}'")
            
def checkInputFiles(anno=None, pangenome=None, fasta=None):
    """
        Checks if the provided input files exist and are of the proper format
    """
    if pangenome is not None:
        if not os.path.exists(pangenome):
            raise FileNotFoundError(f"No such file or directory: '{pangenome}'")

    if anno is not None:
        if not os.path.exists(anno):
            raise FileNotFoundError(f"No such file or directory: '{anno}'")
        checkTsvSanity(anno)

    if fasta is not None:
        if not os.path.exists(fasta):
            raise FileNotFoundError(f"No such file or directory: '{fasta}'")
        checkTsvSanity(fasta)

def cmdLine():

    #need to manually write the description so that it's displayed into groups of subcommands ....
    desc = "\n"
    desc += "  Basic:\n"
    desc += "    workflow      Easy workflow to run a pangenome analysis in one go\n"
    desc += "    panrgp        Easy workflow to run a pangenome analysis with genomic islands and spots of insertion detection\n"
    desc += "  \n"
    desc += "  Expert:\n"
    desc += "    annotate      Annotate genomes\n"
    desc += "    cluster       Cluster proteins in protein families\n"
    desc += "    graph         Create the pangenome graph\n"
    desc += "    partition     Partition the pangenome graph\n"
    desc += "    rarefaction     Compute the rarefaction curve of the pangenome\n"
    desc += "  \n"
    desc += "  Output:\n"
    desc += "    draw          Draw figures representing the pangenome through different aspects\n"
    desc += "    write         Writes 'flat' files representing the pangenome that can be used with other softwares\n"
    desc += "    info          Prints information about a given pangenome graph file\n"
    desc += "  \n"
    desc += "  Regions of genomic Plasticity:\n"
    desc += "    align        aligns a genome or a set of proteins to the pangenome gene families representatives and predict informations from it\n"
    desc += "    rgp          predicts Regions of Genomic Plasticity in the genomes of your pangenome\n"
    desc += "    spot         predicts spots in your pangenome\n"

    parser = argparse.ArgumentParser(description = "Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + pkg_resources.get_distribution("ppanggolin").version)
    subparsers = parser.add_subparsers( metavar = "", dest="subcommand", title="subcommands", description = desc)
    subparsers.required = True#because python3 sent subcommands to hell apparently

    subs = []#subparsers
    subs.append(ppanggolin.annotate.syntaSubparser(subparsers))
    subs.append(ppanggolin.cluster.clusterSubparser(subparsers))
    subs.append(ppanggolin.graph.graphSubparser(subparsers))
    subs.append(ppanggolin.nem.partition.partitionSubparser(subparsers))
    subs.append(ppanggolin.nem.rarefaction.rarefactionSubparser(subparsers))
    subs.append(ppanggolin.workflow.workflow.workflowSubparser(subparsers))
    subs.append(ppanggolin.workflow.panRGP.panRGPSubparser(subparsers))
    subs.append(ppanggolin.figures.figureSubparser(subparsers))
    subs.append(ppanggolin.formats.writeFlat.writeFlatSubparser(subparsers))
    subs.append(ppanggolin.align.alignSubparser(subparsers))
    subs.append(ppanggolin.RGP.genomicIsland.rgpSubparser(subparsers))
    subs.append(ppanggolin.RGP.hotspot.hotspotSubparser(subparsers))
    ppanggolin.info.infoSubparser(subparsers)#not adding to subs because the 'common' options are not needed for this.

    for sub in subs:#add options common to all subcommands
        common = sub._action_groups.pop(1)#get the 'optional arguments' action group.
        common.title = "Common arguments"
        common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(), help = "directory for storing temporary files")
        common.add_argument("--verbose",required=False, type=int,default=1,choices=[0,1,2], help = "Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
        common.add_argument("-c","--cpu",required = False, default = 1,type=int, help = "Number of available cpus")
        common.add_argument('-f', '--force', action="store_true", help="Force writing in output directory and in pangenome output file.")
        sub._action_groups.append(common)
        if (len(sys.argv) == 2 and sub.prog.split()[1] == sys.argv[1]):
            sub.print_help()
            exit(1)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    if args.subcommand == "annotate":
        if args.fasta is None and args.anno is None:
            raise Exception( "You must provide at least a file with the --fasta option to annotate from sequences, or a file with the --gff option to load annotations from.")
    return args

def main():
    args = cmdLine()

    if hasattr(args, "pangenome"):
        checkInputFiles(pangenome = args.pangenome)
    if hasattr(args, "fasta"):
        checkInputFiles(fasta = args.fasta)
    if hasattr(args,"anno"):
        checkInputFiles(anno = args.anno)

    if hasattr(args, "verbose"):
        if args.verbose == 2:
            level = logging.DEBUG#info, debug, warnings and errors
        elif args.verbose == 1:
            level = logging.INFO#info, warnings and errors
        elif args.verbose == 0:
            level = logging.WARNING#only warnings and errors
        logging.basicConfig(stream=sys.stdout, level = level, format = '%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        logging.getLogger().info("Command: "+" ".join([arg for arg in sys.argv]))
        logging.getLogger().info("PPanGGOLiN version: "+pkg_resources.get_distribution("ppanggolin").version)

    if args.subcommand == "annotate":
        ppanggolin.annotate.launch(args)
    elif args.subcommand == "cluster":
        ppanggolin.cluster.launch(args)
    elif args.subcommand == "graph":
        ppanggolin.graph.launch(args)
    elif args.subcommand == "partition":
        ppanggolin.nem.partition.launch(args)
    elif args.subcommand == "workflow":
        ppanggolin.workflow.workflow.launch(args)
    elif args.subcommand == "rarefaction":
        ppanggolin.nem.rarefaction.launch(args)
    elif args.subcommand == "draw":
        ppanggolin.figures.launch(args)
    elif args.subcommand == "write":
        ppanggolin.formats.launch(args)
    elif args.subcommand == "info":
        ppanggolin.info.launch(args)
    elif args.subcommand == "align":
        ppanggolin.align.launch(args)
    elif args.subcommand == "rgp":
        ppanggolin.RGP.genomicIsland.launch(args)
    elif args.subcommand == "spot":
        ppanggolin.RGP.hotspot.launch(args)
    elif args.subcommand == "panrgp":
        ppanggolin.workflow.panRGP.launch(args)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
#coding:utf-8
# PYTHON_ARGCOMPLETE_OK

#default libraries
import sys
if sys.version_info < (3, 6):#minimum is python3.6
    raise AssertionError("Minimum python version to run PPanGGOLiN is 3.6. Your current python version is " + ".".join(map(str,sys.version_info)))
import argparse
import logging
import resource
import pkg_resources
import tempfile
#libraries to be installed
import psutil

try:
    import argcomplete
except ImportError:
    pass

#local modules
import ppanggolin.pangenome
import ppanggolin.nem
import ppanggolin.graph
import ppanggolin.annotate
import ppanggolin.cluster
import ppanggolin.workflow
import ppanggolin.figures
import ppanggolin.formats
import ppanggolin.info
import ppanggolin.align
import ppanggolin.RGP

def requirements():
    """
        Checks if the tools and libraries required for each submodule are installed.
    """
    pass

def cmdLine():

    #need to manually write the description so that it's displayed into groups of subcommands ....
    desc = "\n"
    desc += "  Basic:\n"
    desc += "    workflow      Easy workflow to run a pangenome analysis in one go without parameter tuning\n"
    desc += "  \n"
    desc += "  Expert:\n"
    desc += "    annotate      Annotate genomes\n"
    desc += "    cluster       Cluster proteins in protein families\n"
    desc += "    graph         Create the pangenome graph\n"
    desc += "    partition     Partition the pangenome graph\n"
    desc += "    evolution     Compute the evolution curve of the pangenome\n"
    desc += "  \n"
    desc += "  Output:\n"
    desc += "    draw          Draw figures representing the pangenome through different aspects\n"
    desc += "    write         Writes 'flat' files representing the pangenome that can be used with other softwares\n"
    desc += "    info          Prints information about a given pangenome graph file\n"
    desc += "  \n"
    desc += "  Specific analysis:\n"
    desc += "    align        aligns proteins to the pangenome gene families representatives\n"
    desc += "    rgp          predicts Regions of Genomic Plasticity in the genomes of your pangenome"

    parser = argparse.ArgumentParser(description = "Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v','--version', action='version', version='%(prog)s ' + pkg_resources.get_distribution("ppanggolin").version)
    subparsers = parser.add_subparsers( metavar = "", dest="subcommand", title="subcommands", description = desc)
    subparsers.required = True#because python3 sent subcommands to hell apparently

    subs = []#subparsers
    subs.append(ppanggolin.annotate.syntaSubparser(subparsers))
    subs.append(ppanggolin.cluster.clusterSubparser(subparsers))
    subs.append(ppanggolin.graph.graphSubparser(subparsers))
    subs.append(ppanggolin.nem.partition.partitionSubparser(subparsers))
    subs.append(ppanggolin.nem.evolution.evolutionSubparser(subparsers))
    subs.append(ppanggolin.workflow.workflowSubparser(subparsers))
    subs.append(ppanggolin.figures.figureSubparser(subparsers))
    subs.append(ppanggolin.formats.writeFlat.writeFlatSubparser(subparsers))
    subs.append(ppanggolin.align.alignSubparser(subparsers))
    subs.append(ppanggolin.RGP.rgpSubparser(subparsers))
    ppanggolin.info.infoSubparser(subparsers)#not adding to subs because the 'common' options are not needed for this.

    for sub in subs:#add options common to all subcommands
        common = sub.add_argument_group(title = "Common options")
        common.add_argument("--tmpdir", required=False, type=str, default=tempfile.gettempdir(), help = "directory for storing temporary files")
        common.add_argument("--verbose",required=False, type=int,default=1,choices=[0,1,2], help = "Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
        common.add_argument("-c","--cpu",required = False, default = 1,type=int, help = "Number of available cpus")
        common.add_argument('-f', '--force', action="store_true", help="Force writing in output directory and in pangenome output file.")
        common.add_argument("-se", "--seed", type = int, default = 42, help="seed used to generate random numbers")
        common.add_argument("--memory", required=False, type=int, default=int(4 * psutil.virtual_memory().total / 5), help="Max amount of allowed RAM in Bytes. Default is 4/5 of the system's RAM. Will work only for graph generation part of the program. The C parts and multiprocessed parts might use more without raising an error")

        if len(sys.argv) == 2 and sub.prog.split()[1] == sys.argv[1]:
            sub.print_help()
            exit(1)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    if "argcomplete" in sys.modules:
        argcomplete.autocomplete(parser)

    args = parser.parse_args()
    if args.subcommand == "annotate":
        if args.fasta is None and args.anno is None:
            raise Exception( "You must provide at least a file with the --fasta option to annotate from sequences, or a file with the --gff option to load annotations from.")
    return args

def main():
    args = cmdLine()
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
    # uncomment to save log in file
    # fhandler = logging.FileHandler(filename = args.output + "/PPanGGOLiN.log",mode = "w")
    # fhandler.setFormatter(logging.Formatter(fmt = "%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s", datefmt='%Y-%m-%d %H:%M:%S'))
    # fhandler.setLevel(level if level != logging.WARNING else logging.INFO)
    # logging.getLogger().addHandler(fhandler)
    if hasattr(args,"memory"):
        rsrc = resource.RLIMIT_DATA
        _, hard = resource.getrlimit(rsrc)
        resource.setrlimit(rsrc, (args.memory, hard))  # limiting allowed RAM


    if args.subcommand == "annotate":
        ppanggolin.annotate.launch(args)
    elif args.subcommand == "cluster":
        ppanggolin.cluster.launch(args)
    elif args.subcommand == "graph":
        ppanggolin.graph.launch(args)
    elif args.subcommand == "partition":
        ppanggolin.nem.partition.launch(args)
    elif args.subcommand == "workflow":
        ppanggolin.workflow.launch(args)
    elif args.subcommand == "evolution":
        ppanggolin.nem.evolution.launch(args)
    elif args.subcommand == "draw":
        ppanggolin.figures.launch(args)
    elif args.subcommand == "write":
        ppanggolin.formats.launch(args)
    elif args.subcommand == "info":
        ppanggolin.info.printInfo(args)
    elif args.subcommand == "align":
        ppanggolin.align.launch(args)
    elif args.subcommand == "rgp":
        ppanggolin.RGP.launch(args)

if __name__ == "__main__":
    main()
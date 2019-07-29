#!/usr/bin/env python3
#coding:utf-8

#default libraries
import argparse
from multiprocessing import Pool
import logging
from pathlib import Path

#installed libraries
from tqdm import tqdm

#local libraries
from ppanggolin.annotate import  annotate_organism
from ppanggolin.pangenome import Pangenome
from ppanggolin.utils import read_compressed_or_not, getCurrentRAM
from ppanggolin.formats import writePangenome, readPangenome

def launchAnnotateOrganism(pack):
    return annotate_organism(*pack)

def launch(args):
    if args.fasta is not None:
        #do stuff
        logging.getLogger().debug(f"RAM before doing anything: {getCurrentRAM()}")
        pangenome = Pangenome()
        logging.getLogger().info(f"Reading {args.fasta} the list of organism files")
        
        arguments = []
        for line in read_compressed_or_not(args.fasta):
            elements = [el.strip() for el in line.split("\t")]
            if len(elements)<=1:
                logging.getLogger().error("No tabulation separator found in organisms file")
                exit(1)
            arguments.append((elements[0], elements[1], elements[2:], args.translation_table, args.kingdom, args.norna, args.tmpdir, args.overlap))
        
        logging.getLogger().info(f"Annotating {len(arguments)} genomes using {args.cpu} cpus...")
        with Pool(processes = args.cpu) as p:
            bar = tqdm(range(len(arguments)), unit = "genome")
            for organism in p.imap_unordered(launchAnnotateOrganism, arguments):
                bar.update()
                pangenome.addOrganism(organism)
        logging.getLogger().info("Done annotating genomes")
        pangenome.status["genomesAnnotated"] = "Yes"#the pangenome is now annotated.
        bar.close()
        logging.getLogger().debug(f"RAM at the end of annotating genomes: {getCurrentRAM()}")
        writePangenome(pangenome, Path(args.output + "/" + args.basename ))
    elif args.gff is not None:
        raise NotImplementedError()

def syntaSubparser(subparser):
    parser = subparser.add_parser("annotate",help = "Annotate genomes")
    optional = parser.add_argument_group(title = "Optional arguments")
    optional.add_argument('--overlap', required=False, action='store_false',default=True, help="Use to not remove genes overlapping with RNA features.")
    optional.add_argument("--norna", required=False, action="store_true", default=False, help="Use to avoid annotating RNA features.")
    optional.add_argument("--kingdom",required = False, type = str.lower, default = "bacteria", choices = ["bacteria","archaea"], help = "Kingdom to which the prokaryota belongs to, to know which models to use for rRNA annotation.")
    optional.add_argument("--translation_table",required=False, default="11", help = "Translation table (genetic code) to use.")
    optional.add_argument("--basename",required = False, default = "pangenome", help = "basename for the output file")
    required = parser.add_argument_group(title = "Required arguments", description = "One of the following arguments is required :")
    required.add_argument('--fasta',  required=False, type=str, help="A tab-separated file listing the organism names, and the fasta filepath of its genomic sequence(s) (the fastas can be compressed). One line per organism.")
    required.add_argument('--gff', required=False, type=str, help="A tab-separated file listing the organism names, and the gff filepath of its annotations (the gffs can be compressed). One line per organism.")
    return parser